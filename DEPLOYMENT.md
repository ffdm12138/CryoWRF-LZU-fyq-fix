# Deployment Guide

How to reproduce the CryoWRF deployment used for the Hexi Corridor (河西走廊) simulations
from this repository. Everything needed is in the tree **except the raw input data**
(ERA5 GRIB files and WPS static geography data), which must be downloaded separately.

## Repository layout

| Path | Contents |
|------|----------|
| `WPS/` | WRF Preprocessing System sources (geogrid / ungrib / metgrid), `namelist.wps` (3 domains, 9/3/1 km), `Vtable.ECMWF` |
| `namelist.input` | WRF runtime config: 3 nested domains, SNOWPACK + MYNN + Thompson, demo period 2022-02-01..10 |
| `io.ini` | SNOWPACK configuration (METEOPATH and SNOWPATH both `./input`) |
| `tslist`, `outhist_variables_ALPS.txt` | time-series stations, output variable list |
| `runtime/` | helper scripts: `make_sno_from_met.py`, `run_ungrib.sh`, `pipeline_metgrid_to_real.sh`, snow-event analysis scripts |
| `run/` | WRF static tables (radiation, aerosol, ishmael/P3 lookup tables) — needed at runtime |

## Data dependencies (NOT in this repo — download separately)

1. **ERA5 meteorological forcing** — Copernicus Climate Data Store
   ([Copernicus License](https://cds.climate.copernicus.eu/api/terms/)):
   0.25°, 3-hourly, region 0N–80N, 40E–150E; `pl` + `sfc` GRIB files named
   `0N80N_40E150E_{pl,sfc}_YYYYMMDD.grib`.
2. **WPS static geography data** (`geog_data_path` in `namelist.wps`):
   download from the WRF users' page (e.g. `wps_geog` with SRTM 3s), or let geogrid
   use its bundled low-resolution set.
3. **Toolchain**: gfortran ≥ 9, MPICH 4.x, NetCDF-C 4.9.2 (+ NetCDF-Fortran), HDF5 1.12.3.
   `runtime/run_ungrib.sh` shows the exact library paths used in the reference deployment.

## Build

```bash
# WPS (fast)
cd WPS
./configure          # select dmpar
./compile            # produces geogrid/ungrib/metgrid

# WRF (CryoWRF; the phys/cryowrf coupler is already vendored — no submodule init needed)
cd ..                # repo root
./configure          # dmpar + gfortran
./compile em_real    # produces main/wrf.exe, main/real.exe
```

## Data pipeline (mirrors the reference deployment)

1. `geogrid` → `wpsout/geogrid/geo_em.d0{1,2,3}.nc`
2. `ungrib` (ERA5 pl+sfc in a single run, `Vtable.ECMWF`) → `wpsout/ungrib/FILE:*`
3. `metgrid` → `wpsout/metgrid/met_em.d0{1,2,3}.*.nc`
4. `real.exe` → `wrfinput_d0*`, `wrfbdy_d01`, `wrflowinp_d0*`
5. `python3 runtime/make_sno_from_met.py` → `input/snpack_*.sno` + `input/COSMO.smet`
   (SNOWPACK initialization; run with `sn_start_from_file = .true.`)
6. `mpirun -np 32 ./wrf.exe`

Adapt the paths inside `namelist.wps`, `namelist.input` and `runtime/*.sh` to your machine.

## Gotchas learned in production

- **`io.ini` must contain exactly one `SNOWPATH = input` line.** A duplicate entry in a
  `[Output]` section silently disabled the SNOWPACK coupling: all snow output was -999
  while wrfout looked normal — extremely misleading.
- When inspecting WRF outputs with NCO/netCDF4 tools:
  `export HDF5_USE_FILE_LOCKING=FALSE`
- `ungrib` needs the MPI/NetCDF/HDF5 libraries on `LD_LIBRARY_PATH`
  (see `runtime/run_ungrib.sh`).
- With `use_adaptive_time_step`, the SNOWPACK clock is advanced by the actual WRF
  timestep; a fixed increment assumed a constant dt and drifted ~26% fast.
- WRF binaries embed RPATH to the build libraries; rebuilding is required when the
  toolchain location changes.

## Licensing summary

- WRF + WPS: NCAR/UCAR — public-domain code with attribution request; see `LICENSE.txt`
  and `WPS/README`.
- CryoWRF: GPLv3 — see `phys/cryowrf/LICENSE.txt`.
- MeteoIO / SNOWPACK: see `phys/cryowrf/external/{meteoio,snowpack}/` (License.txt/COPYING).
- ERA5 data and static geography data are **not** redistributed here; obtain them under
  their respective terms (Copernicus / NCAR).
