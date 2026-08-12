#!/usr/bin/env python3
"""
分析 2021年11月 强降雪事件（wrfout 时间序列 + 空间分布）。

用法:
  /home/liguang/miniconda3/bin/python3 runtime/analyze_snow_event.py
  /home/liguang/miniconda3/bin/python3 runtime/analyze_snow_event.py -d 3    # 只分析 d03
  /home/liguang/miniconda3/bin/python3 runtime/analyze_snow_event.py --max   # 只看域最大
"""

import os

# 必须在 import netCDF4 之前设置
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import sys
import glob
import argparse
from datetime import datetime

import numpy as np
import netCDF4 as nc

WRFOUT_DIR = "/wrfrun/2021-2023_hexizoulang/wrfout"


def read_series(domain, var):
    """读取指定域所有 wrfout 文件的变量，返回 (times, data 3D/4D)"""
    files = sorted(glob.glob(os.path.join(WRFOUT_DIR, f"wrfout_d{domain:02d}_*")))
    files = [f for f in files if os.path.isfile(f)]
    if not files:
        return None, None

    times = []
    chunks = []
    for fp in files:
        with nc.Dataset(fp, 'r') as f:
            if var not in f.variables:
                print(f"  [WARN] {var} not in {fp}")
                continue
            tstr = [b''.join(x).decode() for x in f.variables['Times'][:]]
            times.extend(tstr)
            chunks.append(f.variables[var][:].data)

    if not chunks:
        return None, None
    data = np.concatenate(chunks, axis=0)
    return times, data


def analyze_domain(domain, landmask_data=None):
    """分析单个域"""
    print(f"\n{'='*80}")
    print(f"  域 d{domain}（{['', '9km', '3km', '1km'][domain]} 分辨率）")
    print(f"{'='*80}")

    # 读取变量
    times, snownc = read_series(domain, 'SNOWNC')   # 累计降雪 mm
    _, snow = read_series(domain, 'SNOW')           # 雪水当量 kg/m2
    _, snowh = read_series(domain, 'SNOWH')         # 雪深 m
    _, sr = read_series(domain, 'SR')               # 固态降水比例
    _, t2 = read_series(domain, 'T2')               # 2m 温度 K

    if times is None:
        print("  无数据！")
        return

    n = len(times)
    print(f"  时次: {n} (从 {times[0]} 至 {times[-1]})")

    # 陆地掩码（取第一帧）
    # 用 LU_INDEX 过滤水域
    with nc.Dataset(os.path.join(WRFOUT_DIR, f"wrfout_d{domain:02d}_{times[0]}"), 'r') as f:
        lu = f.variables['LU_INDEX'][0].data
    land = (lu != 17) & (lu != 21)  # 非水体
    n_land = int(land.sum())

    # 日累计降雪（SNOWNC 差分）
    print(f"\n  {'时间':<20} {'域均SNOW':>9} {'域均SNOWH':>10} {'最大SNOWH':>10} {'日降雪(max)':>12} {'日降雪(mean)':>12} {'SR>0.8占比':>10}")
    print(f"  {'':20} {'(mm)':>9} {'(cm)':>10} {'(cm)':>10} {'(mm/day)':>12} {'(mm/day)':>12} {'(%)':>10}")
    print(f"  {'-'*83}")

    daily_rows = []
    prev_snownc = None
    prev_time = None

    for i in range(n):
        t = times[i]
        d = datetime.strptime(t, "%Y-%m-%d_%H:%M:%S")

        # 每天 18:00 算日累计（或最后时次）
        if d.hour == 18 or i == n - 1:
            cur_snownc = snownc[i, land].sum() / n_land  # 域均累计
            cur_snow = snow[i, land]
            cur_snowh = snowh[i, land]
            cur_sr = sr[i, land]

            day_snow_max = np.nan
            day_snow_mean = np.nan
            if prev_snownc is not None:
                day_snow_max = np.max(snownc[i, land] - prev_snownc) if i > 0 else np.nan
                # 域均日降雪
                day_snow_mean = cur_snownc - prev_snownc

            sr_frac = float(np.mean(cur_sr > 0.8)) * 100 if cur_sr.size else 0

            date_str = d.strftime("%Y-%m-%d")
            row = {
                'date': date_str,
                'mean_snow': float(np.mean(cur_snow)),
                'mean_snowh_cm': float(np.mean(cur_snowh)) * 100,
                'max_snowh_cm': float(np.max(cur_snowh)) * 100,
                'day_snow_max': day_snow_max if day_snow_max is not None else np.nan,
                'day_snow_mean': day_snow_mean if day_snow_mean is not None else np.nan,
                'sr_frac': sr_frac,
            }
            daily_rows.append(row)
            print(f"  {date_str:<20} {row['mean_snow']:>9.1f} {row['mean_snowh_cm']:>10.2f} "
                  f"{row['max_snowh_cm']:>10.2f} "
                  f"{row['day_snow_max'] if not np.isnan(row['day_snow_max']) else '-':>12}"
                  f"{row['day_snow_mean'] if not np.isnan(row['day_snow_mean']) else '-':>12}"
                  f"{row['sr_frac']:>10.0f}")

            prev_snownc = snownc[i, land].sum() / n_land

    # 找暴雪中心（一次性读经纬度）
    with nc.Dataset(os.path.join(WRFOUT_DIR, f"wrfout_d{domain:02d}_{times[0]}"), 'r') as f:
        lons = f.variables['XLONG'][0].data
        lats = f.variables['XLAT'][0].data

    print(f"\n  🔍 降雪中心定位（SNOWNC 差分最大格点, 3h 步长）")
    for i in range(1, n):
        d = datetime.strptime(times[i], "%Y-%m-%d_%H:%M:%S")
        delta = snownc[i].copy() - snownc[i-1]
        delta[~land] = 0
        idx = np.unravel_index(np.argmax(delta), delta.shape)
        maxval = delta[idx]
        if maxval > 0.5:
            lon = lons[idx[0], idx[1]]
            lat = lats[idx[0], idx[1]]
            print(f"    {times[i]}: 3h降雪 {maxval:6.2f} mm @ ({lat:.3f}°N, {lon:.3f}°E)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--domain', type=int, default=1, help='域 (1/2/3)')
    args = parser.parse_args()

    for dom in ([args.domain] if args.domain in (1, 2, 3) else [1, 2, 3]):
        analyze_domain(dom)


if __name__ == "__main__":
    main()
