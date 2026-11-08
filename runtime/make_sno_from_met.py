#!/usr/bin/env python3
"""REMOTE: canonical simplified CryoWRF cold-profile initialization.

Adapted from the CryoWRF example at commit 73ff27306b42825e6a91640ea41c039ed7ed2712
(test/em_real/input/make_sno_for_land.py; Varun Sharma), and GMD 16, 719 (2023),
section 4.3.2. ERA5/WPS four-layer adapter, not the example's COSMO soil units.
Only writes a NEW profile directory. Does not run WRF/real, modify met fields,
smooth horizontally, or manufacture a meteorological forcing time series.
The case's existing METEOPATH/COSMO.smet is managed separately from SNOWPATH.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from netCDF4 import Dataset, chartostring

ALGORITHM = "cryowrf_official_simple_met_v1"
SNOW_LAYERS = 30
MIN_SNOW_DEPTH = 0.03
ICE_DENSITY = 917.0  # Actual bound SNOWPACK 3.4.5 ElementData::updDensity.
WATER_DENSITY = 1000.0
SOIL_SOLID = 0.7
# SNOWPACK DataClasses.cc:2377-2384 reserves pore space for freezing expansion.
SOIL_WATER_MAX = (1.0 - SOIL_SOLID) * ICE_DENSITY / WATER_DENSITY
SOIL_THICKNESS = (1.89, 0.72, 0.21, 0.07)  # bottom to top, metres
SOIL_NE = (3, 3, 3, 2)
SOIL_TEMP = ("ST100289", "ST028100", "ST007028", "ST000007")
SOIL_MOISTURE = ("SM100289", "SM028100", "SM007028", "SM000007")
FIELDS = ("timestamp Layer_Thick T Vol_Frac_I Vol_Frac_W Vol_Frac_V "
          "Vol_Frac_S Rho_S Conduc_S HeatCapac_S rg rb dd sp mk mass_hoar ne CDot metamo")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def nml_integers(text, key, count):
    # Only the exact integer array is accepted; comments cannot supply values.
    clean = re.sub(r"!.*", "", text)
    m = re.search(r"(?im)^\s*" + re.escape(key) + r"\s*=\s*([\d,\s]+)", clean)
    if not m:
        raise ValueError(f"Missing or unsupported namelist integer array: {key}")
    values = [int(v) for v in re.findall(r"\d+", m.group(1))]
    if len(values) < count:
        raise ValueError(f"{key}: expected {count} explicitly specified entries")
    return values[:count]


def read_domain(path, stamp, expected_shape):
    with Dataset(path) as ds:
        times = np.atleast_1d(chartostring(ds["Times"][:])).tolist()
        if times != [stamp]:
            raise ValueError(f"{path}: Times {times!r} != {[stamp]!r}")
        landclass = str(getattr(ds, "MMINLU", "")).strip()
        if landclass not in {"MODIFIED_IGBP_MODIS_NOAH", "MODIS"}:
            raise ValueError(f"Unsupported land classification {landclass!r}; add a reviewed adapter")
        water_codes = {int(ds.ISWATER), int(ds.ISLAKE)}
        if water_codes != {17, 21}:
            raise ValueError(f"Unexpected MODIS water/lake codes: {water_codes}")
        units = {"SNOW": {"kg m-2"}, "SNOWH": {"m"}, "SKINTEMP": {"K"}}
        units.update({v: {"K"} for v in SOIL_TEMP})
        units.update({v: {"m3 m-3", "m3/m3"} for v in SOIL_MOISTURE})
        names = ("XLAT_M", "XLONG_M", "HGT_M", "SKINTEMP", "LU_INDEX", "SNOW", "SNOWH") + SOIL_TEMP + SOIL_MOISTURE
        data = {}
        for name in names:
            v = ds[name]
            if name in units and str(getattr(v, "units", "")).strip() not in units[name]:
                raise ValueError(f"{name}: unsupported units {getattr(v, 'units', None)!r}")
            a = np.asarray(np.ma.filled(v[0], np.nan), dtype=np.float64)
            if a.shape != expected_shape or not np.isfinite(a).all():
                raise ValueError(f"{name}: nonfinite, masked, or unexpected shape {a.shape}")
            data[name] = a
        if not np.equal(data["LU_INDEX"], np.rint(data["LU_INDEX"])).all():
            raise ValueError("Nonintegral land-use categories")
        for name in ("SKINTEMP",) + SOIL_TEMP:
            if ((data[name] < 100) | (data[name] > 400)).any():
                raise ValueError(f"{name}: implausible absolute temperature")
    active = data["SNOWH"] > MIN_SNOW_DEPTH
    density = np.zeros(expected_shape)
    density[active] = data["SNOW"][active] / data["SNOWH"][active]
    if ((density[active] <= 0) | (density[active] >= ICE_DENSITY)).any():
        raise ValueError("Active snow has invalid SWE/depth density; do not silently clip it")
    data["active"] = active
    data["density"] = density
    data["water"] = np.isin(data["LU_INDEX"], sorted(water_codes))
    return data


def profile_text(domain, y, x, d, stamp):
    snowh = float(d["SNOWH"][y, x]) if d["active"][y, x] else 0.0
    count = SNOW_LAYERS if snowh else 0
    station = f"snpack_{domain}_{x + 1}_{y + 1}"
    water = bool(d["water"][y, x])
    # Preserve the example's surrogate water-bed treatment, with MODIS codes.
    solid = 0.99 if water else SOIL_SOLID
    rho_s, conductivity, heatcap = (1000., .6, 4217.) if water else (2400., 2., 2000.)
    grain, albedo, z0 = (10000., .01, .0002) if water else (1000., .2, .02)
    header = {
        "station_id": station, "station_name": "WRF:met_em",
        "latitude": f"{d['XLAT_M'][y, x]:.12f}",
        "longitude": f"{d['XLONG_M'][y, x]:.12f}",
        "altitude": f"{d['HGT_M'][y, x]:.12f}", "nodata": -999,
        "source": ALGORITHM, "ProfileDate": stamp.replace("_", "T"),
        "HS_Last": f"{snowh:.12f}", "SlopeAngle": 0, "SlopeAzi": 0,
        "nSoilLayerData": 4, "nSnowLayerData": count,
        "SoilAlbedo": albedo, "BareSoil_z0": z0,
        "CanopyHeight": 0, "CanopyLeafAreaIndex": 0,
        "CanopyDirectThroughfall": 1, "WindScalingFactor": 1,
        "ErosionLevel": count, "TimeCountDeltaHS": 0, "fields": FIELDS,
    }
    lines = ["SMET 1.1 ASCII", "[HEADER]"]
    lines.extend(f"{k} = {v}" for k, v in header.items())
    lines.append("[DATA]")
    ts = stamp.replace("_", "T")

    def row(thick, temp, ice, liquid, air, soil, rho, cond, cap, rg, rb, ne):
        phases = [ice, liquid, air, soil]
        if min(phases) < -1e-14 or max(phases) > 1 or abs(sum(phases) - 1) > 1e-12:
            raise ValueError(f"Invalid phases in {station}")
        vals = [thick, temp, ice, liquid, air, soil, rho, cond, cap, rg, rb, 0, 1]
        return ts + " " + " ".join(f"{v:.12f}" for v in vals) + f" 1 0 {ne} 0 0"

    for k, (tv, wv) in enumerate(zip(SOIL_TEMP, SOIL_MOISTURE)):
        temp = max(float(d["SKINTEMP"][y, x]), 278.15) if water else float(d[tv][y, x])
        liquid = 0. if water else float(np.clip(d[wv][y, x], 0, SOIL_WATER_MAX))
        lines.append(row(SOIL_THICKNESS[k], temp, 0, liquid, 1 - solid - liquid,
                         solid, rho_s, conductivity, heatcap, grain, 0, SOIL_NE[k]))
    if count:
        temperatures = np.minimum(np.linspace(d[SOIL_TEMP[-1]][y, x], d["SKINTEMP"][y, x], SNOW_LAYERS), 273.15)
        ice = float(d["density"][y, x]) / ICE_DENSITY
        for k, temp in enumerate(temperatures):
            lines.append(row(snowh / SNOW_LAYERS, temp, ice, 0, 1 - ice, 0,
                             0, 0, 0, .4, .1, 2 if k >= SNOW_LAYERS - 2 else 1))
    return station + ".sno", "\n".join(lines) + "\n"


def domain_statistics(d):
    active, land = d["active"], ~d["water"]
    moisture = np.stack([d[v] for v in SOIL_MOISTURE])[:, land]
    clipped = np.clip(moisture, 0, SOIL_WATER_MAX)
    snow, hs = d["SNOW"], d["SNOWH"]
    return {
        "shape": list(hs.shape), "profiles": int(hs.size),
        "snow_profiles": int(active.sum()), "water_profiles": int(d["water"].sum()),
        "snow_on_water_profiles": int((active & d["water"]).sum()),
        "negative_depth_points": int((hs < 0).sum()),
        "negative_swe_points": int((snow < 0).sum()),
        "positive_depth_removed_by_3cm_cutoff": int(((hs > 0) & ~active).sum()),
        "sum_gridpoint_swe_retained_kg_m2": float(snow[active].sum()),
        "sum_gridpoint_swe_removed_kg_m2": float(snow[(hs > 0) & ~active].sum()),
        "max_snow_depth_m": float(hs[active].max()) if active.any() else 0,
        "soil_land_layer_values": int(moisture.size),
        "soil_moisture_min": float(moisture.min()), "soil_moisture_max": float(moisture.max()),
        "soil_negative_clipped": int((moisture < 0).sum()),
        "soil_above_freezable_capacity_clipped": int((moisture > SOIL_WATER_MAX).sum()),
        "sum_gridpoint_soil_water_change_kg_m2": float(np.sum((clipped - moisture) * np.array(SOIL_THICKNESS)[:, None]) * WATER_DENSITY),
        "subfreezing_wet_soil_rows": int(((np.stack([d[v] for v in SOIL_TEMP])[:, land] < 273.15) & (clipped > 0)).sum()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--date", required=True, help="Exact UTC YYYY-MM-DD_HH:MM:SS, matching namelist start")
    parser.add_argument("--met-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory inside case; never overwrite or skip old files")
    parser.add_argument("-n", "--nworkers", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.nworkers <= 32:
        parser.error("nworkers must be 1..32")
    case, met, out = args.case_dir.resolve(), args.met_dir.resolve(), args.output_dir.resolve()
    if out == case or case not in out.parents or out.exists():
        parser.error("output must be a new child directory of the selected case")
    stamp = datetime.strptime(args.date, "%Y-%m-%d_%H:%M:%S").strftime("%Y-%m-%d_%H:%M:%S")
    nml_path = case / "namelist.input"
    nml = nml_path.read_text()
    ndom = nml_integers(nml, "max_dom", 1)[0]
    starts = [nml_integers(nml, k, ndom) for k in ("start_year", "start_month", "start_day", "start_hour", "start_minute", "start_second")]
    for d in range(ndom):
        if datetime(*(a[d] for a in starts)).strftime("%Y-%m-%d_%H:%M:%S") != stamp:
            parser.error(f"domain {d + 1}: requested date differs from namelist")
    nx, ny = nml_integers(nml, "e_we", ndom), nml_integers(nml, "e_sn", ndom)
    datasets = []
    for domain in range(1, ndom + 1):
        path = met / f"met_em.d{domain:02d}.{stamp}.nc"
        data = read_domain(path, stamp, (ny[domain - 1] - 1, nx[domain - 1] - 1))
        datasets.append((domain, path, data))
    staging = out.with_name(out.name + ".building")
    staging.mkdir(parents=False, exist_ok=False)
    manifest = {
        "algorithm": ALGORITHM, "created_utc": datetime.now(timezone.utc).isoformat(),
        "generator_sha256": sha256(__file__), "namelist_sha256": sha256(nml_path),
        "case": str(case), "output": str(out), "profile_date_utc": stamp,
        "upstream_commit": "73ff27306b42825e6a91640ea41c039ed7ed2712",
        "snow_layers": SNOW_LAYERS, "snow_cutoff_m": MIN_SNOW_DEPTH,
        "density_ice_kg_m3": ICE_DENSITY, "soil_water_max": SOIL_WATER_MAX,
        "precision_decimal_places": 12, "horizontal_smoothing": False,
        "meteo_stub": "Not generated. Preserve/configure METEOPATH separately; SNOWPATH selects this directory.",
        "soil_initial_phase": "Liquid water from volumetric met moisture, ice zero; BUCKET phase change equilibrates later",
        "domains": {}, "files": {},
    }
    for domain, path, data in datasets:
        info = domain_statistics(data)
        info.update({"met_file": str(path), "met_sha256": sha256(path)})
        manifest["domains"][f"d{domain:02d}"] = info
        def write_point(point):
            y, x = point
            name, text = profile_text(domain, y, x, data, stamp)
            raw = text.encode("ascii")
            with (staging / name).open("xb") as f:
                f.write(raw)
            return name, {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
        points = np.ndindex(data["SNOWH"].shape)
        with ThreadPoolExecutor(max_workers=args.nworkers) as pool:
            for name, receipt in pool.map(write_point, points):
                manifest["files"][name] = receipt
        print(json.dumps({"domain": domain, **info}), flush=True)
    # No incomplete directory is published as the requested output.
    manifest["profile_count"] = len(manifest["files"])
    (staging / "initialization-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    staging.rename(out)
    print(json.dumps({"status": "GENERATED_PENDING_INDEPENDENT_VALIDATION", "profiles": manifest["profile_count"], "output": str(out)}), flush=True)


if __name__ == "__main__":
    main()
