#!/usr/bin/env bash
set -euo pipefail

# --- WPS ungrib wrapper for 2021-2023_hexizoulang ---
# ERA5 data via Vtable.ECMWF (unified pressure-level + surface, single pass)
# GRIB data: /wrfdata/metdata/0N80N_40E150E_{pl,sfc}_YYYYMMDD.grib

export LD_LIBRARY_PATH="/home/liguang/CRYOWRF-V1.2.4/libraries_new/mpich/lib:/home/liguang/CRYOWRF-V1.2.4/libraries_new/netcdf/lib:/home/liguang/CRYOWRF-V1.2.4/libraries_new/hdf5/lib:/home/liguang/CRYOWRF-V1/libraries/lib:${LD_LIBRARY_PATH:-}"

cd /wrfrun/2021-2023_hexizoulang/WPS

echo "=== ungrib.exe start: $(date) ==="
./ungrib.exe
echo "=== ungrib.exe exit: $?  $(date) ==="
