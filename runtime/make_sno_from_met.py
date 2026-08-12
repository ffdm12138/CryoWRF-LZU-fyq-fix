#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 met_em 文件生成 SNOWPACK 初始积雪剖面文件（.sno）和 COSMO.smet。
仅使用 met_em 数据源（不依赖 wrfinput）。

功能：
  - 自动从 namelist.input / namelist.wps 读取路径配置
  - 支持多域并行处理
  - 使用共享内存避免大型数组在进程间复制
  - 默认 64 核并行

用法：
  python make_sno_from_met.py                  # 处理所有活跃域，64 核并行
  python make_sno_from_met.py -n 32            # 指定 32 核
  python make_sno_from_met.py -d 1 2           # 仅处理 d01 和 d02
  python make_sno_from_met.py --date 2022-04-25 # 指定日期
  python make_sno_from_met.py --no-cosmo       # 不生成 COSMO.smet
  python make_sno_from_met.py --overwrite      # 覆盖已有文件

@authors: Varun, Guang (原始 make_sno_for_land.py)
@date: 2026-06-04
"""

import os
import sys
import time
import argparse
import multiprocessing as mp
from datetime import datetime, date

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Namelist 解析
# ---------------------------------------------------------------------------

def parse_namelist(filepath):
    """解析 Fortran namelist 文件，返回 {key: value_string} 字典。"""
    result = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!') or line.startswith('#'):
                continue
            if '=' in line:
                key, val = line.split('=', 1)
                result[key.strip()] = val.strip().rstrip(',').strip().strip("'\"")
    return result


def parse_namelist_list(filepath, key):
    """从 namelist 中读取指定 key 的列表值（支持多行逗号分隔）。"""
    vals = []
    found = False
    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('!') or stripped.startswith('#'):
                continue
            if key in stripped and '=' in stripped:
                found = True
                after_eq = stripped.split('=', 1)[1]
                vals.extend([v.strip().strip("'\"") for v in after_eq.split(',')])
                if any(v.endswith('&') or v == '&' for v in vals):
                    vals = [v.rstrip('&').strip() for v in vals if v != '&']
                    continue
                break
            elif found:
                vals.extend([v.strip().strip("'\"") for v in stripped.split(',')])
                if any(v.endswith('&') or v == '&' for v in vals):
                    vals = [v.rstrip('&').strip() for v in vals if v != '&']
                    continue
                break
    return [v for v in vals if v]


def get_domains(config):
    """从 namelist.input 解析活跃域信息。"""
    base = os.path.dirname(os.path.abspath(__file__))
    namelist_input = config.get('namelist_input',
                                os.path.join(base, '..', 'namelist.input'))
    nml = parse_namelist(namelist_input)
    max_dom = int(nml.get('max_dom', '1'))

    def parse_int_list(key, default=None):
        raw = parse_namelist_list(namelist_input, key)
        if raw:
            return [int(x) for x in raw[:max_dom]]
        return default or [1] * max_dom

    return {
        'max_dom': max_dom,
        'e_we': parse_int_list('e_we', [133]),
        'e_sn': parse_int_list('e_sn', [103]),
        'start_year': parse_int_list('start_year', [2022]),
        'start_month': parse_int_list('start_month', [1]),
        'start_day': parse_int_list('start_day', [1]),
        'start_hour': parse_int_list('start_hour', [0]),
    }


def get_met_em_path(config):
    """获取 met_em 文件所在目录路径。"""
    base = os.path.dirname(os.path.abspath(__file__))
    namelist_input = config.get('namelist_input',
                                os.path.join(base, '..', 'namelist.input'))
    nml = parse_namelist(namelist_input)
    run_dir = os.path.dirname(namelist_input)

    # 优先检查常见路径
    candidates = [
        os.path.join(run_dir, 'wpsout', 'metgrid'),
        os.path.join(run_dir, 'met_em'),
        os.path.join(run_dir, 'wps', 'metgrid'),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return c + '/'

    # 从 input_inname 推断
    template = nml.get('input_inname', '')
    if template:
        path = os.path.dirname(template)
        if not os.path.isabs(path):
            path = os.path.join(run_dir, path)
        met_dir = os.path.join(os.path.dirname(path), 'met_em')
        if os.path.isdir(met_dir):
            return met_dir + '/'

    # 从 namelist.wps 读取
    wps_dir = os.environ.get('WPS_DIR', os.path.expanduser('~/CRYOWRF-V1/wps'))
    namelist_wps = os.path.join(wps_dir, 'namelist.wps')
    if not os.path.exists(namelist_wps):
        namelist_wps = os.path.join(base, '..', 'namelist.wps')
    if os.path.exists(namelist_wps):
        nml_wps = parse_namelist(namelist_wps)
        path = nml_wps.get('opt_output_from_metgrid_path', '')
        if path and os.path.isdir(path):
            return path.rstrip('/') + '/'

    return './'


# ---------------------------------------------------------------------------
# 数据读取
# ---------------------------------------------------------------------------

def find_met_em_file(met_dir, domain, dt_obj):
    """查找 met_em 文件，支持多种命名格式。"""
    d = dt_obj.strftime('%Y-%m-%d_%H:%M:%S')
    d_nosec = dt_obj.strftime('%Y-%m-%d_%H:%M')
    patterns = [
        f'met_em.d0{domain}.{d}.nc',
        f'met_em.d0{domain}.{d_nosec}.nc',
        f'met_em.d{domain:02d}.{d}.nc',
    ]
    for pat in patterns:
        fpath = os.path.join(met_dir, pat)
        if os.path.exists(fpath):
            return fpath
    # glob 回退
    import glob
    prefix = os.path.join(met_dir, f'met_em.d0{domain}.')
    matches = sorted(glob.glob(prefix + '*.nc'))
    if matches:
        return matches[0]
    return None


def read_met_em_data(met_dir, domain, dt_obj):
    """
    从 met_em 文件读取所有需要的变量，返回字典。
    坐标/地形/地表温度从 met_em 获取（XLONG_M, XLAT_M, HGT_M, SKINTEMP）。
    """
    fpath = find_met_em_file(met_dir, domain, dt_obj)
    if fpath is None:
        raise FileNotFoundError(
            f"找不到 met_em 文件: d0{domain}, {dt_obj.isoformat()}")

    ds = xr.open_dataset(fpath)

    lon = ds['XLONG_M'].values[0, :, :]
    lat = ds['XLAT_M'].values[0, :, :]
    hgt = ds['HGT_M'].values[0, :, :]
    skin_temp = ds['SKINTEMP'].values[0, :, :]

    t_soil = np.concatenate([
        ds['ST100289'].values,
        ds['ST028100'].values,
        ds['ST007028'].values,
        ds['ST000007'].values
    ], axis=0)

    w_soil = np.concatenate([
        ds['SM100289'].values,
        ds['SM028100'].values,
        ds['SM007028'].values,
        ds['SM000007'].values
    ], axis=0)

    lu_index = ds['LU_INDEX'].values[0, :, :]
    snow = ds['SNOW'].values[0, :, :]

    # ERA5 met_em 没有 SNOWH（雪深），从 SWE 估算：密度 ~250 kg/m³
    if 'SNOWH' in ds.variables:
        snowh = ds['SNOWH'].values[0, :, :]
    else:
        rho_snow_default = 250.0  # kg/m³, 中等密度雪
        snowh = snow / rho_snow_default  # SWE (kg/m²) → 深度 (m)

    snowh[snowh <= 0.01] = 0.0
    t_soil_top = ds['ST000007'].values[0, :, :]

    ds.close()
    return {
        'lon': lon, 'lat': lat, 'hgt': hgt, 'skin_temp': skin_temp,
        't_soil': t_soil, 'w_soil': w_soil, 'lu_index': lu_index,
        'snow': snow, 'snowh': snowh, 't_soil_top': t_soil_top,
    }


# ---------------------------------------------------------------------------
# .sno 文件写入（单点）
# ---------------------------------------------------------------------------

def write_sno_file(domain, loc_lat, loc_lon, in_i, in_j, loc_hgt,
                   loc_snowh, soil_thickness, soil_T, soil_w,
                   soil_v, soil_s, soil_ne, t_snow,
                   vol_i, vol_v, no_snow_lays,
                   soil_dens, hcond, hcap, grain_r, alb, z0,
                   dt_obj, soil_i):
    """为单个网格点生成 .sno 文件。"""
    fn_out = f"./input/snpack_{domain}_{in_j + 1}_{in_i + 1}.sno"
    n_snow_lays = 10 if loc_snowh > 0 else 0

    with open(fn_out, 'w') as f:
        f.write("SMET 1.1 ASCII\n")
        f.write("[HEADER]\n")
        f.write(f"station_id       = snpack_{domain}_{in_j + 1}_{in_i + 1}\n")
        f.write("station_name     = WRF:met_em\n")
        f.write(f"latitude         = {loc_lat:.4f}\n")
        f.write(f"longitude        = {loc_lon:.4f}\n")
        f.write(f"altitude         = {loc_hgt:.2f}\n")
        f.write("nodata           = -999\n")
        f.write(f"source           = met_em_to_sno; make_sno_from_met.py, "
                f"{datetime.now().strftime('%Y-%m-%d')}\n")
        f.write(f"ProfileDate      = {dt_obj.isoformat()[:19]}\n")
        f.write(f"HS_Last          = {loc_snowh:.6f}\n")
        f.write("SlopeAngle       = 0.00\n")
        f.write("SlopeAzi         = 0.00\n")
        f.write("nSoilLayerData   = 4\n")
        f.write(f"nSnowLayerData   = {n_snow_lays}\n")
        f.write(f"SoilAlbedo       = {alb}\n")
        f.write(f"BareSoil_z0      = {z0}\n")
        f.write("CanopyHeight     = 0.00\n")
        f.write("CanopyLeafAreaIndex = 0.00\n")
        f.write("CanopyDirectThroughfall = 1.00\n")
        f.write("WindScalingFactor = 1.00\n")
        f.write(f"ErosionLevel     = {n_snow_lays}\n")
        f.write("TimeCountDeltaHS = 0.000000\n")
        f.write("fields           = timestamp Layer_Thick  T  Vol_Frac_I  "
                "Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S "
                "HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo\n")
        f.write("[DATA]\n")

        ts = dt_obj.isoformat()[:19]

        # 4 层土壤
        for m in range(4):
            f.write(
                f"{ts}   "
                f"{soil_thickness[m]:6.3f}   "
                f"{soil_T[m]:7.3f}   "
                f"{soil_i[m]:5.3f}   "
                f"{soil_w[m]:5.3f}   "
                f"{soil_v[m]:5.3f}   "
                f"{soil_s[m]:5.3f}   "
                f"{soil_dens:8.3f}   "
                f"{hcond:6.3f}   "
                f"{hcap:8.3f}   "
                f"{grain_r:9.3f}   "
                f"{0.0:9.3f}   "
                f"{0.0:5.3f}   "
                f"{1.0:5.3f}   "
                f"1   "
                f"0.0   "
                f"{soil_ne[m]}   "
                f"0.0   "
                f"0.0\n"
            )

        # 10 层雪
        if loc_snowh > 0:
            for m in range(10):
                f.write(
                    f"{ts}   "
                    f"{loc_snowh / 10.0:6.3f}   "
                    f"{t_snow[m]:7.3f}   "
                    f"{vol_i:5.3f}   "
                    f"{0.0:5.3f}   "
                    f"{vol_v:5.3f}   "
                    f"{0.0:5.3f}   "
                    f"{0.0:8.3f}   "
                    f"{0.0:6.3f}   "
                    f"{0.0:8.3f}   "
                    f"{0.4:9.3f}   "
                    f"{0.1:9.3f}   "
                    f"{0.0:5.3f}   "
                    f"{1.0:5.3f}   "
                    f"1   "
                    f"0.0   "
                    f"{int(no_snow_lays[m])}   "
                    f"0.0   "
                    f"0.0\n"
                )

    return fn_out


# ---------------------------------------------------------------------------
# 共享内存辅助
# ---------------------------------------------------------------------------

def make_shared(shape, dtype='float32'):
    """创建共享内存数组并返回 numpy 视图。"""
    n = int(np.prod(shape))
    shm = mp.Array('f', n)
    arr = np.frombuffer(shm.get_obj(), dtype='float32').reshape(shape)
    arr[:] = 0
    return shm, arr


def init_worker(shared_arrays, grid_info, dt_str, overwrite):
    """Pool 初始化函数，将共享数据绑定到全局变量。"""
    global _data, _grid, _dt, _overwrite
    _data = shared_arrays
    _grid = grid_info
    _dt = datetime.fromisoformat(dt_str)
    _overwrite = overwrite


def process_grid_point(args):
    """处理单个网格点（Pool worker 调用）。"""
    i, j = args

    # 检查是否已存在
    domain = _grid['domain']
    fn = f"./input/snpack_{domain}_{j + 1}_{i + 1}.sno"
    if not _overwrite and os.path.exists(fn):
        return None

    loc_lat = _data['lat'][i, j]
    loc_lon = _data['lon'][i, j]
    loc_hgt = _data['hgt'][i, j]
    loc_snowh = _data['snowh'][i, j]
    soil_T = _data['t_soil'][:, i, j].copy()
    t_snow = _data['t_snow'][:, i, j].copy()

    # 土壤体积分数（与 make_sno_for_land.py 一致：干土初始条件）
    soil_s = np.full(4, 0.7)
    soil_v = np.full(4, 0.3)
    soil_w = np.zeros(4)
    soil_i = np.zeros(4)

    # 雪密度 → 体积分数（与原始 make_sno_for_land.py 一致）
    if loc_snowh > 0:
        rho_snow = _data['snow'][i, j] / loc_snowh
        vol_i = rho_snow / 918.0
        vol_v = 1.0 - vol_i
    else:
        vol_i, vol_v = 0.0, 1.0

    # 默认土壤/地表参数
    soil_dens, hcond, hcap = 2400.0, 2.0, 2000.0
    grain_r, alb, z0 = 1000.0, 0.2, 0.02

    # 水面点特殊处理
    lu = _data['lu_index'][i, j]
    if lu == 17 or lu == 21:
        soil_dens, hcond, hcap = 1000.0, 0.6, 4217.0
        grain_r, alb, z0 = 10000.0, 0.01, 0.0002
        soil_v[:] = 0.01
        soil_w[:] = 0.0
        soil_s[:] = 0.99
        soil_i[:] = 0.0
        soil_T[:] = np.maximum(_data['skin_temp'][i, j], 278.15)

    return write_sno_file(
        domain, loc_lat, loc_lon, i, j, loc_hgt, loc_snowh,
        _grid['soil_layer_thickness'], soil_T, soil_w, soil_v, soil_s,
        _grid['soil_layer_ne'], t_snow, vol_i, vol_v, _grid['no_snow_lays'],
        soil_dens, hcond, hcap, grain_r, alb, z0, _dt, soil_i)


# ---------------------------------------------------------------------------
# COSMO.smet 生成
# ---------------------------------------------------------------------------

def generate_cosmo_smet(domain, dt_obj, met_dir=None, timestep=10, n_hours=2):
    """
    从 met_em 数据生成 COSMO.smet。

    直接从 met_em 文件提取气象变量（TA, RH, VW, DW, ISWR, ILWR, PSUM），
    生成 SNOWPACK 初始化时需要的 COSMO.smet 文件。
    station_id 与 .sno 文件匹配（snpack_<domain>_1_1）。

    Parameters
    ----------
    domain : int
        域编号
    dt_obj : datetime
        模拟开始时间
    met_dir : str, optional
        met_em 文件目录（默认自动查找）
    timestep : int
        数据时间步长（秒），默认 10
    n_hours : int
        生成数据时长（小时），默认 2
    """
    import glob as _glob

    if met_dir is None:
        config = {}
        met_dir = get_met_em_path(config)

    # 查找 met_em 文件
    fpath = find_met_em_file(met_dir, domain, dt_obj)
    if fpath is None:
        print(f"[WARN] 找不到 met_em 文件 (d0{domain}, {dt_obj})，跳过 COSMO.smet 生成")
        return

    ds = xr.open_dataset(fpath)

    # 提取第一个网格点 (0,0) 的气象变量
    t2 = float(ds['TT'].values[0, 0, 0, 0])        # 温度 (K)
    rh = float(ds['RH'].values[0, 0, 0, 0])         # 相对湿度 (%)
    u10 = float(ds['UU'].values[0, 0, 0, 0])        # U 风速 (m/s)
    v10 = float(ds['VV'].values[0, 0, 0, 0])        # V 风速 (m/s)

    # 辐射（可能不在 met_em 中）
    swdown = 0.0
    lwdown = 0.0
    if 'SWDOWN' in ds.variables:
        swdown = float(ds['SWDOWN'].values[0, 0, 0])
    if 'GLW' in ds.variables:
        lwdown = float(ds['GLW'].values[0, 0, 0])

    lat = float(ds['XLAT_M'].values[0, 0, 0])
    lon = float(ds['XLONG_M'].values[0, 0, 0])
    alt = float(ds['HGT_M'].values[0, 0, 0])
    ds.close()

    # 风速和风向
    vw = float(np.sqrt(u10**2 + v10**2))
    dw = float(np.degrees(np.arctan2(-u10, -v10)) % 360)

    # RH: met_em 存储为百分比 (0-100)，SNOWPACK 需要小数 (0-1)
    rh_frac = rh / 100.0 if rh > 1.0 else rh

    iswr = max(swdown, 0.0)
    ilwr = max(lwdown, 0.0)

    station_id = f"snpack_{domain}_1_1"

    # 写入 COSMO.smet
    n_steps = int(n_hours * 3600 / timestep)
    with open('./input/COSMO.smet', 'w') as f:
        f.write("SMET 1.1 ASCII\n")
        f.write("[HEADER]\n")
        f.write(f"station_id       = {station_id}\n")
        f.write(f"station_name     = WRF_SNPACK_d0{domain}\n")
        f.write(f"latitude         = {lat:.6f}\n")
        f.write(f"longitude        = {lon:.6f}\n")
        f.write(f"altitude         = {alt:.1f}\n")
        f.write("nodata           = -999\n")
        f.write("tz               = 8\n")
        f.write("fields           = timestamp TA RH VW DW ISWR ILWR PSUM\n")
        f.write("[DATA]\n")
        for i in range(n_steps):
            t = dt_obj + __import__('datetime').timedelta(seconds=i * timestep)
            ts = t.strftime("%Y-%m-%dT%H:%M:%S")
            f.write(f"{ts}  {t2:.2f}  {rh_frac:.4f}  {vw:.1f}  {dw:.0f}  "
                    f"{iswr:.0f}  {ilwr:.0f}  0.0000\n")

    print(f"  COSMO.smet 已生成: station_id={station_id}, "
          f"TA={t2:.1f}K, RH={rh_frac:.3f}, VW={vw:.1f}m/s, "
          f"日期={dt_obj.isoformat()[:19]}, {n_steps} 步")


# ---------------------------------------------------------------------------
# 主处理逻辑
# ---------------------------------------------------------------------------

def write_sno_files(met_dir, domain, dt_obj, n_workers=64, overwrite=False,
                    gen_cosmo=False):
    """
    为指定域从 met_em 数据生成所有 .sno 文件。

    Parameters
    ----------
    met_dir : str
        met_em 文件所在目录
    domain : int
        域编号（1, 2, ...）
    dt_obj : datetime
        生成初始场的日期时间
    n_workers : int
        并行进程数（默认 64）
    overwrite : bool
        是否覆盖已有的 .sno 文件
    gen_cosmo : bool
        是否生成 COSMO.smet
    """
    print(f"读取 met_em 数据: d0{domain}, {dt_obj.isoformat()}")
    data = read_met_em_data(met_dir, domain, dt_obj)
    print(f"  网格尺寸: {data['lon'].shape}")

    # 预计算雪相关派生量
    snowh = data['snowh']
    snow = data['snow']
    t_soil_top = data['t_soil_top']
    skin_temp = data['skin_temp']

    with np.errstate(divide='ignore', invalid='ignore'):
        rho_snow = np.where(snowh > 0, snow / snowh, 0.0)
    vol_i = rho_snow / 918.0
    vol_v = 1.0 - vol_i

    # 30 层雪温度插值（向量化）
    t_snow = np.linspace(t_soil_top, skin_temp, num=30, axis=0)
    t_snow = np.minimum(t_snow, 273.15)

    # 雪层有限元数（10 层）
    no_snow_lays = np.ones(10)
    no_snow_lays[-2:] = 2

    soil_layer_thickness = np.array([1.89, 0.72, 0.21, 0.07])
    soil_layer_ne = np.array([3, 3, 3, 2])

    nx, ny = data['lon'].shape

    # --- 共享内存（避免进程间复制大型数组） ---
    print("初始化共享内存...")
    shm_lon, arr_lon = make_shared((nx, ny))
    shm_lat, arr_lat = make_shared((nx, ny))
    shm_hgt, arr_hgt = make_shared((nx, ny))
    shm_skin, arr_skin = make_shared((nx, ny))
    shm_ts, arr_ts = make_shared((4, nx, ny))
    shm_ws, arr_ws = make_shared((4, nx, ny))
    shm_lu, arr_lu = make_shared((nx, ny), 'int32')
    shm_snow, arr_snow = make_shared((nx, ny))
    shm_snowh, arr_snowh = make_shared((nx, ny))
    shm_tsnow, arr_tsnow = make_shared((30, nx, ny))

    arr_lon[:] = data['lon']
    arr_lat[:] = data['lat']
    arr_hgt[:] = data['hgt']
    arr_skin[:] = skin_temp
    arr_ts[:] = data['t_soil']
    arr_ws[:] = data['w_soil']
    arr_lu[:] = data['lu_index']
    arr_snow[:] = snow
    arr_snowh[:] = snowh
    arr_tsnow[:] = t_snow

    shared_arrays = {
        'lon': arr_lon, 'lat': arr_lat, 'hgt': arr_hgt,
        'skin_temp': arr_skin, 't_soil': arr_ts, 'w_soil': arr_ws,
        'lu_index': arr_lu, 'snow': arr_snow, 'snowh': arr_snowh,
        't_snow': arr_tsnow,
    }

    grid_info = {
        'domain': domain,
        'nx': nx, 'ny': ny,
        'soil_layer_thickness': soil_layer_thickness,
        'soil_layer_ne': soil_layer_ne,
        'no_snow_lays': no_snow_lays,
    }

    # --- 构造网格点列表 ---
    if not overwrite:
        grid_points = [
            (i, j)
            for i in range(nx) for j in range(ny)
            if not os.path.exists(f"./input/snpack_{domain}_{j + 1}_{i + 1}.sno")
        ]
        skip_count = nx * ny - len(grid_points)
        if skip_count > 0:
            print(f"  跳过 {skip_count} 个已有 .sno 文件")
    else:
        grid_points = [(i, j) for i in range(nx) for j in range(ny)]

    total = len(grid_points)
    print(f"  需要生成 {total} 个 .sno 文件（并行 {n_workers} 核）")

    if total == 0:
        print("  所有 .sno 文件已存在，无需生成")
        if gen_cosmo:
            generate_cosmo_smet(domain, dt_obj, met_dir=met_dir)
        return

    # --- 并行处理 ---
    t0 = time.time()
    with mp.Pool(
        processes=n_workers,
        initializer=init_worker,
        initargs=(shared_arrays, grid_info, dt_obj.isoformat(), overwrite)
    ) as pool:
        results = pool.map(process_grid_point, grid_points)

    created = sum(1 for r in results if r is not None)
    elapsed = time.time() - t0
    print(f"  完成: 生成 {created} 个 .sno 文件, 耗时 {elapsed:.1f}s")

    # --- COSMO.smet ---
    if gen_cosmo:
        generate_cosmo_smet(domain, dt_obj, met_dir=met_dir)


# ---------------------------------------------------------------------------
# 命令行入口
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='从 met_em 文件生成 SNOWPACK 初始积雪剖面（.sno）')
    parser.add_argument('-n', '--nworkers', type=int, default=64,
                        help='并行进程数（默认: 64）')
    parser.add_argument('-d', '--domains', type=int, nargs='+', default=None,
                        help='处理指定域（默认: 所有活跃域）')
    parser.add_argument('--date', type=str, default=None,
                        help='指定日期（格式: YYYY-MM-DD，默认: namelist 中的起始日期）')
    parser.add_argument('--no-cosmo', action='store_true',
                        help='不生成 COSMO.smet')
    parser.add_argument('--overwrite', action='store_true',
                        help='覆盖已有的 .sno 文件')
    args = parser.parse_args()

    config = {}
    domains_info = get_domains(config)
    met_dir = get_met_em_path(config)
    max_dom = domains_info['max_dom']

    dom_list = args.domains if args.domains else list(range(1, max_dom + 1))

    if args.date:
        dt_obj = datetime.strptime(args.date, '%Y-%m-%d')
    else:
        dt_obj = datetime(
            domains_info['start_year'][0],
            domains_info['start_month'][0],
            domains_info['start_day'][0],
            domains_info['start_hour'][0],
        )

    # 切换到工作目录根（河西走廊: 脚本在 runtime/，上一层即为工作目录）
    work_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(work_dir)
    print(f"工作目录: {work_dir}")
    print(f"met_em 目录: {met_dir}")
    print(f"日期: {dt_obj.isoformat()}")
    print(f"处理域: {dom_list}")
    print(f"并行核数: {args.nworkers}")
    print()

    for domain in dom_list:
        print(f"{'=' * 60}")
        print(f"处理域 d0{domain}")
        print(f"{'=' * 60}")
        t0 = time.time()
        write_sno_files(
            met_dir, domain, dt_obj,
            n_workers=args.nworkers,
            overwrite=args.overwrite,
            gen_cosmo=(not args.no_cosmo),
        )
        print(f"域 d0{domain} 完成, 总耗时 {time.time() - t0:.1f}s\n")

    print("全部完成。")


if __name__ == '__main__':
    main()
