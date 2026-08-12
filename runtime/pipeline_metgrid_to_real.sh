#!/bin/bash
# 河西走廊 metgrid → sno → real 自动化流水线
# 运行: bash /wrfrun/2021-2023_hexizoulang/runtime/pipeline_metgrid_to_real.sh

set -e
WORK=/wrfrun/2021-2023_hexizoulang
LOG=$WORK/logs/pipeline.log
PYTHON=/home/liguang/miniconda3/bin/python3

exec >> $LOG 2>&1
echo "=========================================="
echo "PIPELINE START: $(date)"
echo "=========================================="

# ---- Step 0: 等待 metgrid 完成 ----
echo "[0/3] 等待 metgrid.exe 完成..."
while pgrep -f "metgrid.exe" > /dev/null 2>&1; do
  d01=$(ls $WORK/wpsout/metgrid/met_em.d01.*.nc 2>/dev/null | wc -l)
  d02=$(ls $WORK/wpsout/metgrid/met_em.d02.*.nc 2>/dev/null | wc -l)
  d03=$(ls $WORK/wpsout/metgrid/met_em.d03.*.nc 2>/dev/null | wc -l)
  total=$((d01+d02+d03))
  echo "[$(date '+%H:%M:%S')] metgrid: d01:$d01 d02:$d02 d03:$d03 = $total/17523 ($((total*100/17523))%)"
  sleep 300
done
echo "metgrid.exe 已退出: $(date)"

# ---- 验证 metgrid 输出 ----
D01=$(ls $WORK/wpsout/metgrid/met_em.d01.*.nc 2>/dev/null | wc -l)
D02=$(ls $WORK/wpsout/metgrid/met_em.d02.*.nc 2>/dev/null | wc -l)
D03=$(ls $WORK/wpsout/metgrid/met_em.d03.*.nc 2>/dev/null | wc -l)
echo "metgrid 最终: d01:$D01 d02:$D02 d03:$D03"
if [ "$D01" -ne 5841 ] || [ "$D02" -ne 5841 ] || [ "$D03" -ne 5841 ]; then
  echo "FATAL: metgrid 文件数不完整！"
  exit 1
fi
grep -q "Successful completion" $WORK/logs/metgrid.stdout.log && echo "metgrid: SUCCESS" || { echo "metgrid: 无成功标记"; }

# ---- Step 1: make_sno_from_met.py ----
echo ""
echo "[1/3] 生成 .sno + COSMO.smet: $(date)"
cd $WORK
$PYTHON runtime/make_sno_from_met.py -n 64
echo "make_sno_from_met: DONE ($(date))"

SNO_COUNT=$(ls $WORK/input/snpack_*.sno 2>/dev/null | wc -l)
echo ".sno 文件数: $SNO_COUNT"
[ -f "$WORK/input/COSMO.smet" ] && echo "COSMO.smet: OK" || echo "COSMO.smet: MISSING"

# ---- Step 2: real.exe ----
echo ""
echo "[2/3] 运行 real.exe: $(date)"
cd $WORK
mpirun -np 8 ./real.exe > $WORK/logs/real.stdout.log 2> $WORK/logs/real.stderr.log
REAL_EXIT=$?
echo $REAL_EXIT > $WORK/logs/real.exitcode

if [ $REAL_EXIT -eq 0 ] && grep -q "SUCCESS COMPLETE" rsl.out.0000 2>/dev/null; then
  echo "real.exe: SUCCESS ($(date))"
  ls -lh $WORK/wrfinput_d0* $WORK/wrfbdy_d01 2>/dev/null
else
  echo "real.exe: FAILED (exit=$REAL_EXIT)"
  tail -20 $WORK/logs/real.stderr.log
  exit 1
fi

# ---- 完成 ----
echo ""
echo "=========================================="
echo "PIPELINE COMPLETE: $(date)"
echo "d01:$D01 d02:$D02 d03:$D03 met_em"
echo "sno: $SNO_COUNT files"
echo "wrfinput: $(ls $WORK/wrfinput_d0* 2>/dev/null | wc -l)"
echo "wrfbdy: $(ls $WORK/wrfbdy_d01 2>/dev/null)"
echo "READY FOR WRF"
echo "=========================================="
