#!/usr/bin/env python3
"""
Interpolate (nearest-neighbour) .sno spinup files to CRYOWRF d01 grid,
with physics-based layer merging (mirrors SNOWPACK Aggregate::aggregate()).

Source:  /capstor/scratch/cscs/rengbers/sno_restart_08052024/VIR_1_<j>_<i>.sno
Target:  CRYOWRF_ANTARCTICA_PEA_Lorena/WRF/wrfinput_d01  (121x121, 27km polar stereo)
Output:  ./snpack_1_<j+1>_<i+1>.sno  (run from CRYOWRF_ANTARCTICA_PEA_Lorena/WRF/)

For each target ISLTYP==8 (glacier ice) cell, finds the nearest source .sno file
by lat/lon, merges layers thinner than MIN_LAYER_THICK using SNOWPACK aggregation
rules, then writes a new file with updated station_id/lat/lon/altitude header.
Prints the great-circle distance for every match and a summary at the end.
"""
import os, glob, re, math
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree
from netCDF4 import Dataset

SNO_DIR         = "/capstor/scratch/cscs/rengbers/sno_restart_08052024"
WRFINPUT        = "/capstor/scratch/cscs/rengbers/CRYOWRF_ANTARCTICA_PEA_Lorena/WRF/wrfinput_d01"
DD_NUM          = 1
MIN_LAYER_THICK = 0.01  # merge layers thinner than 1 cm (SNOWPACK aggregation logic)
MAX_SNOW_LAYERS = 100   # secondary cap: keep at most this many layers after merging
OUT_DIR         = "."  # output directory (run from WRF/)


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two lat/lon points."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat/2)**2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon/2)**2)
    return 2 * R * math.asin(math.sqrt(a))


# --- 1. Parse source .sno headers ---
print("Scanning source .sno files ...")
src_lats, src_lons, src_paths = [], [], []
for f in sorted(glob.glob(os.path.join(SNO_DIR, "VIR_1_*.sno"))):
    lat = lon = None
    with open(f) as fh:
        for line in fh:
            if line.startswith("latitude"):
                lat = float(line.split("=")[1])
            elif line.startswith("longitude"):
                lon = float(line.split("=")[1])
            elif line.startswith("[DATA]"):
                break
    if lat is not None and lon is not None:
        src_lats.append(lat)
        src_lons.append(lon)
        src_paths.append(f)

src_lats = np.array(src_lats)
src_lons = np.array(src_lons)
print(f"  {len(src_paths)} source files loaded.")
tree = KDTree(np.column_stack((src_lons, src_lats)))


# ---------------------------------------------------------------------------
# Snow-layer merging — mirrors SNOWPACK Aggregate::aggregate() second-run pass
# and DataClasses::mergeElements() property combination formulas.
#
# .sno column layout (0-based):
#   0  timestamp (string, Julian date "YYYY-MM-DDTHH:MM:SS")
#   1  Layer_Thick   → sum
#   2  T             → thickness-weighted average
#   3  Vol_Frac_I    → thickness-weighted average
#   4  Vol_Frac_W    → thickness-weighted average
#   5  Vol_Frac_V    → recalculate: 1 - θ_ice - θ_water - θ_soil
#   6  Vol_Frac_S    → thickness-weighted average
#   7  Rho_S         → thickness-weighted average
#   8  Conduc_S      → thickness-weighted average
#   9  HeatCapac_S   → thickness-weighted average
#  10  rg            → ice-content-weighted average
#  11  rb            → ice-content-weighted average
#  12  dd            → ice-content-weighted average
#  13  sp            → ice-content-weighted average
#  14  mk            → lower layer value; add upper century-digit if upper>=100 and lower<100
#  15  mass_hoar     → upper layer value
#  16  ne            → always 1
#  17  CDot          → ice-content-weighted average
#  18  metamo        → thickness-weighted average
# ---------------------------------------------------------------------------

# Thresholds translated from SNOWPACK Aggregate.cc / ElementData.h
_LIMIT_DRY = 0.001   # θ_water below this → "dry" layer
_DIFF_JUL  = 1.0     # max age difference (days) for dry-layer merging


def _julian_days(ts_str):
    """Convert .sno timestamp string to ordinal float (days)."""
    try:
        dt = datetime.strptime(ts_str.strip(), "%Y-%m-%dT%H:%M:%S")
        return dt.toordinal() + dt.hour / 24.0 + dt.minute / 1440.0
    except ValueError:
        return 0.0


def _can_merge(lower, upper):
    """
    Translated from SNOWPACK mergeThinLayer() (ElementData version).
    Returns True if two adjacent layers are compatible for merging.
    lower = deeper (older) layer; upper = shallower (newer) layer.
    Row index 0 in snow_rows is the deepest layer.
    """
    theta_w_lo = lower[4]
    theta_w_up = upper[4]
    mk_lo = int(lower[14]) % 100
    mk_up = int(upper[14]) % 100

    # Dry layers: check water content difference, age, and grain class
    if theta_w_lo < _LIMIT_DRY or theta_w_up < _LIMIT_DRY:
        if abs(theta_w_lo - theta_w_up) > _LIMIT_DRY:
            return False
        if abs(_julian_days(lower[0]) - _julian_days(upper[0])) > 2 * _DIFF_JUL:
            return False
        if mk_lo != mk_up:
            return False
    else:
        # Both wet: check grain size and sphericity thresholds
        rg_lo = lower[10]; sp_lo = lower[13]
        rg_up = upper[10]; sp_up = upper[13]
        if rg_lo < 0.75 or sp_lo < 0.5 or rg_up < 0.75 or sp_up < 0.5:
            if mk_lo != mk_up:
                return False
        else:
            if not (9 < mk_lo < 13 or 19 < mk_lo < 23):
                if mk_lo != mk_up:
                    return False
    return True


def _merge_two(lower, upper):
    """
    Translated from SNOWPACK DataClasses::mergeElements().
    lower = deeper layer, upper = shallower layer.
    Returns merged layer (list with same layout).
    """
    L_lo  = lower[1]
    L_up  = upper[1]
    L_new = L_lo + L_up

    # Ice content (kg/m²-equivalent: θ_ice × thickness)
    ice_lo  = lower[3] * L_lo
    ice_up  = upper[3] * L_up
    ice_sum = ice_lo + ice_up

    # Start from lower layer (inherits position in snowpack)
    merged = list(lower)

    # Thickness
    merged[1] = L_new

    # Thickness-weighted averages
    for col in [2, 3, 4, 6, 7, 8, 9, 18]:
        merged[col] = (lower[col] * L_lo + upper[col] * L_up) / L_new

    # Ice-content-weighted averages (grain properties)
    if ice_sum > 0:
        for col in [10, 11, 12, 13, 17]:
            merged[col] = (lower[col] * ice_lo + upper[col] * ice_up) / ice_sum

    # Recalculate air volume fraction
    merged[5] = max(0.0, 1.0 - merged[3] - merged[4] - merged[6])

    # mk: keep lower's grain class; add upper's century-digit if present
    mk_lo = int(lower[14])
    mk_up = int(upper[14])
    if mk_up >= 100 and mk_lo < 100:
        merged[14] = (mk_up // 100) * 100 + mk_lo
    else:
        merged[14] = mk_lo

    # Timestamp and mass_hoar from upper (newer/surface) layer
    merged[0]  = upper[0]
    merged[15] = upper[15]

    # ne always 1 after merge
    merged[16] = 1.0

    return merged


def merge_thin_layers(text, min_thick=MIN_LAYER_THICK, max_layers=MAX_SNOW_LAYERS):
    """
    Replace thin snow layers with physics-based merged layers, then cap at
    max_layers as a secondary safeguard.

    Mirrors SNOWPACK Aggregate::aggregate() second-run (hard-floor) pass:
    - iterate bottom-to-top
    - if a layer is thin (< min_thick) OR sticky flag is set, attempt merge
    - if layers are compatible OR force is True, merge and re-check same index
    - otherwise set sticky flag to propagate to next pair
    - after each full pass restart if anything changed
    - finally handle any thin top layer
    """
    m_soil = int(re.search(r"(?m)^nSoilLayerData\s*=\s*(\d+)", text).group(1))
    m_snow = int(re.search(r"(?m)^nSnowLayerData\s*=\s*(\d+)", text).group(1))

    data_pos = text.find("[DATA]")
    eol      = text.find("\n", data_pos)
    header   = text[:eol + 1]
    raw_rows = [l for l in text[eol + 1:].splitlines() if l.strip()]

    soil_rows = raw_rows[:m_soil]
    snow_raw  = raw_rows[m_soil:m_soil + m_snow]

    # Parse snow rows: col 0 = timestamp string, cols 1..N = floats
    def parse_row(line):
        parts = line.split()
        return [parts[0]] + [float(x) for x in parts[1:]]

    snow = [parse_row(r) for r in snow_raw]

    # --- Merge pass (bottom-to-top, sticky flag) ---
    flag    = False
    changed = True
    while changed:
        changed = False
        flag    = False
        i = 0
        while i < len(snow) - 1:
            lower = snow[i]
            upper = snow[i + 1]
            force = lower[1] < min_thick
            if force or flag:
                if _can_merge(lower, upper) or force:
                    snow[i] = _merge_two(lower, upper)
                    snow.pop(i + 1)
                    flag    = False
                    changed = True
                    continue    # re-examine same index with new partner
                else:
                    flag = True
            else:
                flag = False
            i += 1

        # Handle thin top layer (merge downward into layer below)
        if len(snow) > 1 and snow[-1][1] < min_thick:
            snow[-2] = _merge_two(snow[-2], snow[-1])
            snow.pop(-1)
            changed = True

    # --- Secondary cap ---
    if len(snow) > max_layers:
        snow = snow[-max_layers:]    # keep topmost (surface) layers

    # --- Serialise rows ---
    def fmt_row(row):
        ts   = row[0]
        nums = row[1:]
        # reproduce original formatting: fixed-width fields
        parts = [ts]
        for v in nums:
            # Use integer formatting for cols that are integer-valued (mk, ne)
            if v == int(v) and abs(v) < 1e6:
                parts.append(f"{int(v):6d}")
            else:
                parts.append(f"{v:12.6f}")
        return " ".join(parts)

    snow_text = "\n".join(fmt_row(r) for r in snow) + "\n"

    # --- Update header fields ---
    n_snow_new = len(snow)
    hs_new     = sum(r[1] for r in snow)

    header = re.sub(r"(?m)^nSnowLayerData\s*=.*$",
                    f"nSnowLayerData   = {n_snow_new}", header)
    header = re.sub(r"(?m)^HS_Last\s*=.*$",
                    f"HS_Last          = {hs_new:.6f}", header)

    m = re.search(r"(?m)^ErosionLevel\s*=\s*(\d+)", header)
    if m:
        old_erosion = int(m.group(1))
        new_erosion = max(0, min(old_erosion, n_snow_new - 1))
        header = re.sub(r"(?m)^ErosionLevel\s*=.*$",
                        f"ErosionLevel     = {new_erosion}", header)

    return header + "\n".join(soil_rows + [fmt_row(r) for r in snow]) + "\n"


# --- 2. Target grid ---
print("Reading wrfinput_d01 ...")
wrf     = Dataset(WRFINPUT, 'r')
lat_wrf = wrf.variables["XLAT"][0,:,:]
lon_wrf = wrf.variables["XLONG"][0,:,:]
isltyp  = wrf.variables["ISLTYP"][0,:,:]
hgt     = wrf.variables["HGT"][0,:,:]
wrf.close()

# --- 3. Interpolate and write ---
max_i, max_j = lat_wrf.shape
distances_km = []
n_written = 0

for i in range(max_i):
    for j in range(max_j):
        if isltyp[i, j] == 8:
            tlat = float(lat_wrf[i, j])
            tlon = float(lon_wrf[i, j])
            thgt = float(hgt[i, j])
            fn_out = f"{OUT_DIR}/snpack_{DD_NUM}_{j+1}_{i+1}.sno"
            new_id = f"snpack_{DD_NUM}_{j+1}_{i+1}"

            _, idx = tree.query([tlon, tlat])
            dist_km = haversine_km(tlat, tlon, src_lats[idx], src_lons[idx])
            distances_km.append(dist_km)

            with open(src_paths[idx]) as fh:
                text = fh.read()

            text = re.sub(r"(?m)^station_id\s*=.*$",
                          f"station_id       = {new_id}", text)
            text = re.sub(r"(?m)^latitude\s*=.*$",
                          f"latitude         = {tlat:.6f}", text)
            text = re.sub(r"(?m)^longitude\s*=.*$",
                          f"longitude        = {tlon:.6f}", text)
            text = re.sub(r"(?m)^altitude\s*=.*$",
                          f"altitude         = {thgt:.2f}", text)
            if "WindScalingFactor" not in text:
                text = re.sub(r"(?m)^(CanopyDirectThroughfall\s*=.*$)",
                              r"\1\nWindScalingFactor = 1.0", text)
            else:
                text = re.sub(r"(?m)^WindScalingFactor\s*=.*$",
                              "WindScalingFactor = 1.0", text)

            text = merge_thin_layers(text)

            with open(fn_out, 'w') as fh:
                fh.write(text)

            print(f"  ({i:3d},{j:3d})  →  {os.path.basename(src_paths[idx])}   "
                  f"dist={dist_km:6.1f} km   →  {fn_out}")

            n_written += 1

# --- 4. Summary ---
print(f"\nDone. {n_written} files written.")
if distances_km:
    print(f"  Mean match distance : {np.mean(distances_km):.1f} km")
    print(f"  Max  match distance : {np.max(distances_km):.1f} km")
    print(f"  Cells > 50 km away  : {sum(d > 50 for d in distances_km)}")
