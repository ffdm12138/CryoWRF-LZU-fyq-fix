#!/usr/bin/env python3
"""
从 ERA5 met_em 文件中识别强降雪时段。

方法：
  追踪 SNOW (积雪水当量, kg/m² ≈ mm SWE) 在相邻 3 小时间隔的变化。
  ΔSNOW > 0 表示积雪增加（降雪 > 融雪）。
  结合近地面气温 (TAVGSFC) 过滤暖期降雨事件。

输出：
  - 日累计降雪量时间序列 (CSV)
  - Top-N 强降雪事件列表
  - 月度降雪总量统计

用法：
  /home/liguang/miniconda3/bin/python3 runtime/find_snow_events.py
  /home/liguang/miniconda3/bin/python3 runtime/find_snow_events.py -n 32   # 并行
"""

import os
import sys
import glob
import argparse
from datetime import datetime, timedelta
from multiprocessing import Pool
from collections import defaultdict

import numpy as np
import netCDF4 as nc

# 配置
METGRID_DIR = "/wrfrun/2021-2023_hexizoulang/wpsout/metgrid"
DOMAIN = "d01"
SNOW_VAR = "SNOW"
TEMP_VAR = "SKINTEMP"  # 地表温度 (K)，随时间变化
LAND_VAR = "LANDMASK"
SNOW_THRESHOLD_KGM2 = 0.1  # 最小 ΔSNOW 变化 (kg/m²)，过滤噪声
TEMP_SNOW_THRESHOLD = 275.15  # 2°C，低于此温度认为降雪可能性大


def get_met_files():
    """获取排序后的 met_em 文件列表"""
    pattern = os.path.join(METGRID_DIR, f"met_em.{DOMAIN}.*.nc")
    files = sorted(glob.glob(pattern))
    return files


def read_snow_and_temp(filepath):
    """读取单个文件的 SNOW, TAVGSFC, LANDMASK 的域平均"""
    try:
        with nc.Dataset(filepath, 'r') as f:
            snow = f.variables[SNOW_VAR][0, :, :].data  # (south_north, west_east)
            temp = f.variables[TEMP_VAR][0, :, :].data
            landmask = f.variables[LAND_VAR][0, :, :].data

        # 只统计陆地格点
        land = (landmask > 0.5) & (~np.isnan(snow)) & (~np.isnan(temp))
        if land.sum() == 0:
            return {
                'mean_snow': np.nan,
                'mean_temp': np.nan,
                'n_land': 0,
                'max_snow': np.nan,
                'snow_positive_pct': np.nan,
            }

        return {
            'mean_snow': float(np.mean(snow[land])),
            'mean_temp': float(np.mean(temp[land])),
            'n_land': int(np.sum(land)),
            'max_snow': float(np.max(snow[land])),
            # 积雪覆盖率 (>1 kg/m²)
            'snow_cover_pct': float(np.sum((snow[land] > 1.0)) / land.sum() * 100),
        }
    except Exception as e:
        print(f"  [ERROR] {os.path.basename(filepath)}: {e}", file=sys.stderr)
        return None


def parse_time_from_filename(filepath):
    """从文件名解析时间"""
    basename = os.path.basename(filepath)
    # met_em.d01.2021-07-01_00:00:00.nc
    time_str = basename.split('.', 2)[-1].replace('.nc', '')
    return datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")


def process_files(files, nproc=16):
    """并行读取所有文件"""
    print(f"读取 {len(files)} 个 met_em 文件 (domain={DOMAIN}, {nproc} 核)...")
    with Pool(nproc) as pool:
        results = pool.map(read_snow_and_temp, files)

    times = [parse_time_from_filename(f) for f in files]

    data = []
    for t, r in zip(times, results):
        if r is not None:
            data.append((t, r))

    data.sort(key=lambda x: x[0])
    print(f"成功读取 {len(data)}/{len(files)} 个文件")
    return data


def compute_daily_snowfall(data):
    """
    计算日累计降雪量。

    方法:
      ΔSNOW(t) = SNOW(t) - SNOW(t-3h)
      如果 ΔSNOW > 0 且 T < 2°C → 计为降雪
      日累计 = Σ(max(ΔSNOW, 0)) over 8 time steps
    """
    daily = defaultdict(lambda: {
        'total_snowfall': 0.0,      # 日累计 ΔSNOW > 0 (kg/m²)
        'cold_snowfall': 0.0,       # 冷期降雪 (T < 2°C)
        'max_hourly_snowfall': 0.0, # 单时次最大
        'n_snow_steps': 0,          # 降雪时次
        'n_cold_steps': 0,
        'mean_temp': 0.0,
        'mean_snow_cover': 0.0,
        'max_snow_cover': 0.0,
        'max_snow': 0.0,
        'n_obs': 0,
    })

    hourly_events = []  # 每 3 小时的事件

    for i in range(1, len(data)):
        t_prev, r_prev = data[i - 1]
        t_curr, r_curr = data[i]

        if np.isnan(r_curr['mean_snow']) or np.isnan(r_prev['mean_snow']):
            continue

        delta_snow = r_curr['mean_snow'] - r_prev['mean_snow']
        temp = r_curr['mean_temp']
        date_key = t_curr.strftime("%Y-%m-%d")

        daily[date_key]['n_obs'] += 1
        daily[date_key]['mean_temp'] += temp
        daily[date_key]['mean_snow_cover'] += r_curr['snow_cover_pct']
        daily[date_key]['max_snow_cover'] = max(
            daily[date_key]['max_snow_cover'], r_curr['snow_cover_pct']
        )
        daily[date_key]['max_snow'] = max(
            daily[date_key]['max_snow'], r_curr['max_snow']
        )

        if delta_snow > SNOW_THRESHOLD_KGM2:
            daily[date_key]['total_snowfall'] += delta_snow
            daily[date_key]['max_hourly_snowfall'] = max(
                daily[date_key]['max_hourly_snowfall'], delta_snow
            )
            daily[date_key]['n_snow_steps'] += 1

            if temp < TEMP_SNOW_THRESHOLD:
                daily[date_key]['cold_snowfall'] += delta_snow
                daily[date_key]['n_cold_steps'] += 1

            hourly_events.append({
                'time': t_curr,
                'delta_snow': delta_snow,
                'temp': temp,
                'snow_cover_pct': r_curr['snow_cover_pct'],
            })

    # 归一化日平均
    for d in daily:
        n = daily[d]['n_obs']
        if n > 0:
            daily[d]['mean_temp'] /= n
            daily[d]['mean_snow_cover'] /= n

    return daily, hourly_events


def find_top_events(daily, hourly_events, top_n=20):
    """识别强降雪事件（连续降雪日）"""
    # 按日降雪量排序
    sorted_days = sorted(daily.items(), key=lambda x: x[1]['total_snowfall'], reverse=True)

    # 找出连续降雪的时段（事件）
    all_dates = sorted(daily.keys())
    events = []
    event_start = None
    event_end = None
    event_snow = 0.0

    for date_str in all_dates:
        snow = daily[date_str]['total_snowfall']
        if snow > 0.5:  # 有效降雪日（>0.5 mm SWE/day）
            if event_start is None:
                event_start = date_str
            event_end = date_str
            event_snow += snow
        else:
            if event_start is not None:
                dur = (datetime.strptime(event_end, "%Y-%m-%d") -
                       datetime.strptime(event_start, "%Y-%m-%d")).days + 1
                events.append({
                    'start': event_start,
                    'end': event_end,
                    'duration_days': dur,
                    'total_snowfall': round(event_snow, 2),
                    'intensity': round(event_snow / dur, 2),
                })
                event_start = None
                event_snow = 0.0

    # 处理最后一个事件
    if event_start is not None:
        dur = (datetime.strptime(event_end, "%Y-%m-%d") -
               datetime.strptime(event_start, "%Y-%m-%d")).days + 1
        events.append({
            'start': event_start,
            'end': event_end,
            'duration_days': dur,
            'total_snowfall': round(event_snow, 2),
            'intensity': round(event_snow / dur, 2),
        })

    events.sort(key=lambda x: x['total_snowfall'], reverse=True)

    # 找出 3 小时强降雪时刻
    top_hourly = sorted(hourly_events, key=lambda x: x['delta_snow'], reverse=True)[:top_n]

    return sorted_days[:top_n], events, top_hourly


def monthly_summary(daily):
    """月度降雪统计"""
    monthly = defaultdict(lambda: {
        'total_snowfall': 0.0,
        'n_snow_days': 0,
        'n_days': 0,
        'mean_temp': 0.0,
        'mean_snow_cover': 0.0,
    })

    for date_str, d in daily.items():
        month_key = date_str[:7]  # YYYY-MM
        monthly[month_key]['total_snowfall'] += d['total_snowfall']
        monthly[month_key]['n_days'] += 1
        monthly[month_key]['mean_temp'] += d['mean_temp']
        monthly[month_key]['mean_snow_cover'] += d['mean_snow_cover']
        if d['n_snow_steps'] > 0:
            monthly[month_key]['n_snow_days'] += 1

    for m in monthly:
        n = monthly[m]['n_days']
        monthly[m]['mean_temp'] = round(monthly[m]['mean_temp'] / n - 273.15, 1)
        monthly[m]['mean_snow_cover'] = round(monthly[m]['mean_snow_cover'] / n, 1)
        monthly[m]['total_snowfall'] = round(monthly[m]['total_snowfall'], 2)

    return monthly


def main():
    parser = argparse.ArgumentParser(description="从 ERA5 met_em 识别强降雪时段")
    parser.add_argument('-n', '--nproc', type=int, default=16, help='并行核数')
    parser.add_argument('-o', '--output', type=str,
                        default='/wrfrun/2021-2023_hexizoulang/runtime/snow_events_output.csv',
                        help='日降雪量输出 CSV')
    parser.add_argument('--top', type=int, default=20, help='显示 Top-N 事件')
    args = parser.parse_args()

    files = get_met_files()
    if not files:
        print("未找到 met_em 文件！", file=sys.stderr)
        sys.exit(1)

    print(f"找到 {len(files)} 个 met_em 文件 ({files[0]} ... {files[-1]})")

    # Step 1: 读取
    data = process_files(files, args.nproc)

    # Step 2: 计算日降雪量
    print("计算日降雪量...")
    daily, hourly_events = compute_daily_snowfall(data)

    # Step 3: 识别强降雪事件
    print("识别强降雪事件...")
    top_days, events, top_hourly = find_top_events(daily, hourly_events, args.top)

    # Step 4: 月度统计
    monthly = monthly_summary(daily)

    # === 输出 ===
    print("\n" + "=" * 90)
    print("                     ERA5 河西走廊 — 强降雪事件识别结果")
    print("=" * 90)

    # 月度统计表
    print(f"\n{'─' * 75}")
    print(f"  📊 月度降雪统计")
    print(f"{'─' * 75}")
    print(f"  {'月份':<10} {'降雪总量':>10} {'降雪日数':>8} {'均温(°C)':>10} {'积雪覆盖率':>10}")
    print(f"  {'':10} {'(mm SWE)':>10} {'(天)':>8} {'':10} {'(%)':>10}")
    print(f"{'─' * 75}")
    for m in sorted(monthly.keys()):
        d = monthly[m]
        print(f"  {m:<10} {d['total_snowfall']:>10.2f} {d['n_snow_days']:>8} "
              f"{d['mean_temp']:>10.1f} {d['mean_snow_cover']:>10.1f}")

    # 强降雪事件
    print(f"\n{'─' * 90}")
    print(f"  🔥 强降雪事件（连续降雪日合并）")
    print(f"{'─' * 90}")
    print(f"  {'排名':<6} {'开始日期':<12} {'结束日期':<12} {'持续(天)':>8} {'累计降雪':>10} {'日均强度':>10}")
    print(f"  {'':6} {'':12} {'':12} {'':>8} {'(mm SWE)':>10} {'(mm/day)':>10}")
    print(f"{'─' * 90}")
    for i, ev in enumerate(events[:args.top], 1):
        print(f"  {i:<6} {ev['start']:<12} {ev['end']:<12} {ev['duration_days']:>8} "
              f"{ev['total_snowfall']:>10.2f} {ev['intensity']:>10.2f}")

    # 最强单日
    print(f"\n{'─' * 75}")
    print(f"  📅 最强单日降雪 (Top-{min(args.top, 15)})")
    print(f"{'─' * 75}")
    print(f"  {'排名':<6} {'日期':<12} {'日降雪量':>10} {'3h最大':>10} {'均温':>10} {'积雪覆盖率':>10}")
    print(f"  {'':6} {'':12} {'(mm SWE)':>10} {'(mm/3h)':>10} {'(°C)':>10} {'(%)':>10}")
    print(f"{'─' * 75}")
    for i, (d, v) in enumerate(top_days[:15], 1):
        temp_c = v['mean_temp'] - 273.15
        print(f"  {i:<6} {d:<12} {v['total_snowfall']:>10.2f} {v['max_hourly_snowfall']:>10.2f} "
              f"{temp_c:>10.1f} {v['max_snow_cover']:>10.1f}")

    # 最强 3 小时
    print(f"\n{'─' * 90}")
    print(f"  ⚡ 最强 3 小时降雪时刻")
    print(f"{'─' * 90}")
    print(f"  {'排名':<6} {'时间':<22} {'3h降雪量':>10} {'温度':>10} {'积雪覆盖率':>10}")
    print(f"  {'':6} {'':22} {'(mm SWE)':>10} {'(°C)':>10} {'(%)':>10}")
    print(f"{'─' * 90}")
    for i, ev in enumerate(top_hourly[:15], 1):
        temp_c = ev['temp'] - 273.15
        print(f"  {i:<6} {ev['time'].strftime('%Y-%m-%d %H:%M'):<22} "
              f"{ev['delta_snow']:>10.3f} {temp_c:>10.1f} {ev['snow_cover_pct']:>10.1f}")

    # 保存 CSV
    with open(args.output, 'w') as f:
        f.write("date,total_snowfall_mmSWE,cold_snowfall_mmSWE,n_snow_steps,max_hourly_mmSWE,"
                "mean_temp_C,max_snow_cover_pct,mean_snow_cover_pct,max_snow_mmSWE\n")
        for date_str in sorted(daily.keys()):
            d = daily[date_str]
            temp_c = d['mean_temp'] - 273.15 if d['mean_temp'] > 0 else np.nan
            f.write(f"{date_str},{d['total_snowfall']:.4f},{d['cold_snowfall']:.4f},"
                    f"{d['n_snow_steps']},{d['max_hourly_snowfall']:.4f},"
                    f"{temp_c:.1f},{d['max_snow_cover']:.2f},{d['mean_snow_cover']:.2f},"
                    f"{d['max_snow']:.2f}\n")

    print(f"\n✅ 日降雪量时间序列已保存至: {args.output}")
    print(f"   共 {len(daily)} 天数据")


if __name__ == "__main__":
    main()
