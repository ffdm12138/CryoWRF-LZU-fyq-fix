#!/usr/bin/env python3
"""
Interpolate (nearest-neighbour) .sno spinup files to CRYOWRF d01 grid.

Source:  /capstor/scratch/cscs/rengbers/sno_restart_08052024/VIR_1_<j>_<i>.sno
Target:  CRYOWRF_ANTARCTICA_PEA_Lorena/WRF/wrfinput_d01  (121x121, 27km polar stereo)
Output:  ./snpack_1_<j+1>_<i+1>.sno  (run from CRYOWRF_ANTARCTICA_PEA_Lorena/WRF/)

For each target ISLTYP==8 (glacier ice) cell, finds the nearest source .sno file
by lat/lon and writes a new file with updated station_id/lat/lon/altitude header.
Prints the great-circle distance for every match and a summary at the end.
"""
import os, glob, re, math
import numpy as np
from scipy.spatial import KDTree
from netCDF4 import Dataset

SNO_DIR         = "/capstor/scratch/cscs/rengbers/sno_restart_08052024"
WRFINPUT        = "/capstor/scratch/cscs/rengbers/CRYOWRF_ANTARCTICA_PEA_Lorena/WRF/wrfinput_d01"
DD_NUM          = 1
MAX_SNOW_LAYERS = 100   # keep only the top N snow layers (surface layers)
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


def truncate_snow_layers(text, max_layers):
    """Keep only the top max_layers snow layers (last rows in DATA section)."""
    m_soil = int(re.search(r"(?m)^nSoilLayerData\s*=\s*(\d+)", text).group(1))
    m_snow = int(re.search(r"(?m)^nSnowLayerData\s*=\s*(\d+)", text).group(1))
    if m_snow <= max_layers:
        return text
    n_removed = m_snow - max_layers
    data_pos = text.find("[DATA]")
    eol = text.find("\n", data_pos)
    header = text[:eol + 1]
    rows = [l for l in text[eol + 1:].splitlines() if l.strip()]
    soil_rows = rows[:m_soil]
    snow_rows = rows[m_soil:][-max_layers:]
    text = header + "\n".join(soil_rows + snow_rows) + "\n"
    text = re.sub(r"(?m)^nSnowLayerData\s*=.*$",
                  f"nSnowLayerData   = {max_layers}", text)
    # Update HS_Last to match the kept layers
    hs = sum(float(r.split()[1]) for r in snow_rows if r.strip())
    text = re.sub(r"(?m)^HS_Last\s*=.*$",
                  f"HS_Last          = {hs:.6f}", text)
    # Update ErosionLevel: 0-based index into the snow layer array
    m = re.search(r"(?m)^ErosionLevel\s*=\s*(\d+)", text)
    if m:
        old_erosion = int(m.group(1))
        new_erosion = max(0, old_erosion - n_removed)
        new_erosion = min(new_erosion, max_layers - 1)
        text = re.sub(r"(?m)^ErosionLevel\s*=.*$",
                      f"ErosionLevel     = {new_erosion}", text)
    return text


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

            text = truncate_snow_layers(text, MAX_SNOW_LAYERS)

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
