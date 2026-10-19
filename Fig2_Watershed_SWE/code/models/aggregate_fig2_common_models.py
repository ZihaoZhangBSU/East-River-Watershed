"""Reproduce fixed-footprint Fig. 2 model means and audit archived summaries.

No model simulations are run. ``aggregate_fixed`` operates on existing daily
Noah-grid spatial states, never on native-grid basin means. The command line
imports the already spatially remapped, independently raw-spot-checked October
2026 extraction and preserves its exact saved means. The audit rejects any
record with incomplete common weighted coverage. Full spatial caches are not
distributed with the publication package, and this import is not a fresh raw
model extraction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_calendar(frame):
    """Require every displayed day exactly once, including leap days."""
    frame = frame.copy()
    frame["date"] = pd.to_datetime(frame["date"], errors="raise")
    if frame["date"].duplicated().any():
        raise ValueError("Duplicate dates")
    if not frame["date"].is_monotonic_increasing:
        raise ValueError("Dates must be increasing")
    if set(frame["wy"].unique()) != set(range(2018, 2027)):
        raise ValueError("Expected WY2018 through WY2026")
    for wy, part in frame.groupby("wy"):
        end = f"{wy}-06-30" if wy == 2026 else f"{wy}-09-30"
        expected = pd.date_range(f"{wy-1}-10-01", end, freq="D")
        if not pd.DatetimeIndex(part["date"]).equals(expected):
            raise ValueError(f"Incomplete or shifted calendar for WY{wy}")
    return frame


def aggregate_fixed(swe_i, swe_n, weights, fixed_mask=None):
    """Area-weight existing Noah-grid fields on one fixed footprint.

    Missing/negative SWE at any retained target cell is an error, not a reason
    to silently renormalize a daily mean. Zeros are valid. Fractional watershed
    areas are the weights; equal-cell weights are never substituted.
    """
    w = np.asarray(weights, dtype=np.float64)
    if w.ndim != 2 or not np.isfinite(w).all() or (w < 0).any():
        raise ValueError("Weights must be a finite, nonnegative 2D array")
    mask = w > 0 if fixed_mask is None else np.asarray(fixed_mask, dtype=bool)
    if mask.shape != w.shape or (mask & (w <= 0)).any() or not mask.any():
        raise ValueError("Invalid fixed footprint")
    outputs = []
    for a in (swe_i, swe_n):
        a = np.asarray(a, dtype=np.float64)
        if a.ndim != 3 or a.shape[1:] != w.shape:
            raise ValueError("Expected daily states with shape (date, row, col)")
        retained = a[:, mask]
        if not np.isfinite(retained).all() or (retained < 0).any():
            raise ValueError("Nonfinite or negative SWE on fixed model footprint")
        outputs.append(np.sum(retained * w[mask], axis=1) / np.sum(w[mask]))
    if len(outputs[0]) != len(outputs[1]):
        raise ValueError("Models have different date counts")
    return tuple(outputs)


def coverage_precision_audit(weights):
    """Check stored float32 coverage can detect even one missing positive cell.

    The original extractor used ``nw[joint_finite].sum()/nw.sum()``. Test
    omissions with the exact dtype and boolean-index reduction order, rather
    than assume a saved 1.0 could not conceal a tiny boundary-cell omission.
    """
    w = np.asarray(weights)
    positive = w > 0
    denominator = w.sum()
    ratios = []
    for j, k in np.argwhere(positive):
        good = positive.copy()
        good[j, k] = False
        ratios.append(float(w[good].sum() / denominator))
    if not ratios or max(ratios) >= 1:
        raise ValueError("Coverage precision cannot resolve one missing cell")
    return {
        "weight_dtype": str(w.dtype),
        "positive_weight_cells": int(positive.sum()),
        "minimum_positive_weight": float(w[positive].min()),
        "float32_sum_fractional_cells": float(denominator),
        "float64_sum_fractional_cells": float(np.sum(w, dtype=np.float64)),
        "full_positive_cell_ratio": float(w[positive].sum()/denominator),
        "maximum_ratio_after_one_positive_cell_omitted": max(ratios),
        "single_cell_omissions_tested": len(ratios),
        "invisible_single_cell_omissions": int(sum(r == 1 for r in ratios)),
    }


def annual_metrics(new, old):
    """Annual displayed and common October--June means, agreement and changes."""
    old = check_calendar(old)
    new = check_calendar(new)
    if not old["date"].equals(new["date"]):
        raise ValueError("Native and common-grid calendars differ")
    rows = []
    for wy in range(2018, 2027):
        for window in ("displayed", "oct_jun"):
            keep = new["wy"].eq(wy)
            if window == "oct_jun":
                keep &= new["date"].dt.month.isin([10,11,12,1,2,3,4,5,6])
            row = {"wy": wy, "window": window, "n_days": int(keep.sum())}
            for label, frame in (("old_native", old), ("new_common_grid", new)):
                i = frame.loc[keep, "swe_i_mm"].to_numpy(float)
                n = frame.loc[keep, "swe_n_mm"].to_numpy(float)
                mi, mn = float(i.mean()), float(n.mean())
                # fsum provides an independent scalar summation check.
                if not math.isclose(mi, math.fsum(i)/len(i), abs_tol=1e-10):
                    raise AssertionError("Independent mean check failed")
                rr=float(np.corrcoef(i,n)[0,1])
                scalar_r=math.fsum((x-mi)*(y-mn) for x,y in zip(i,n))/math.sqrt(math.fsum((x-mi)**2 for x in i)*math.fsum((y-mn)**2 for y in n))
                if not math.isclose(rr,scalar_r,abs_tol=1e-12):
                    raise AssertionError("Independent scalar correlation check failed")
                row.update({
                    f"{label}_mean_i_mm": mi,
                    f"{label}_mean_n_mm": mn,
                    f"{label}_noah_minus_isnobal_percent": 100*(mn/mi-1),
                    f"{label}_r": rr,
                    f"{label}_mb_n_minus_i_mm": float(np.mean(n-i)),
                })
            for tag in ("i", "n"):
                row[f"mean_{tag}_change_mm"] = row[f"new_common_grid_mean_{tag}_mm"]-row[f"old_native_mean_{tag}_mm"]
            rows.append(row)
    return pd.DataFrame(rows)


def import_verified(args):
    src = Path(args.source_dir)
    dst = Path(args.output_dir)
    dst.mkdir(parents=True, exist_ok=True)
    core = check_calendar(pd.read_csv(src/"core_watershed_daily.csv"))
    archived_coverage = core["valid_common_area_fraction"].to_numpy()
    if not (archived_coverage == archived_coverage[0]).all() or archived_coverage.min() < 1-2*np.finfo(np.float32).eps:
        raise ValueError("Archived coverage changes or is below float32 rounding tolerance")
    values = core[["swe_i_mm", "swe_n_mm"]].to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Invalid archived means")
    validation = pd.read_csv(src/"cache_raw_validation.csv")
    if len(validation) != 18 or set(validation["wy"]) != set(range(2018,2027)):
        raise ValueError("Expected two archived raw checks per year")
    if (validation[["isnobal_swe_max_difference_mm", "noah_swe_max_difference_mm"]].abs().to_numpy() > .002).any():
        raise ValueError("Archived cache/raw checks failed")
    methods = json.loads((src/"core_methods.json").read_text())
    if methods["grid_shape"] != [66,54] or "D23" not in methods["state_support"]:
        raise ValueError("Wrong archived extraction")
    with np.load(src/"static_noah.npz", allow_pickle=False) as z:
        static = {k:z[k] for k in z.files}
    w = static["weights"]
    precision = coverage_precision_audit(w)
    grid_meta = json.loads(Path(args.grid_metadata).read_text())["Noah-MP"]
    a,b,c,d,e,f = grid_meta["transform"][:6]
    if not np.allclose(static["x"], c+a*(np.arange(len(static["x"]))+.5), atol=1e-6):
        raise ValueError("Grid x centers differ")
    if not np.allclose(static["y"], f+e*(np.arange(len(static["y"]))+.5), atol=1e-6):
        raise ValueError("Grid y centers differ")
    np.savez_compressed(dst/"grid_common.npz", **static, fixed_model_mask=w>0, fixed_domain_mask=w>0,
                        transform=np.asarray(grid_meta["transform"]),
                        affine=np.asarray(grid_meta["transform"]), crs=np.asarray(grid_meta["crs"]))
    selected = core[["date","wy","swe_i_mm","swe_n_mm"]].copy()
    selected["fixed_model_area_fraction"] = 1.0
    selected["daily_common_valid_area_fraction"] = core["valid_common_area_fraction"]
    selected["fixed_model_cells"] = int((w>0).sum())
    selected["state_hour_utc"] = 23
    selected.to_csv(dst/"Fig2_common_grid_daily.csv",index=False,date_format="%Y-%m-%d",float_format="%.17g")
    native_path = Path(args.native_daily)
    native = pd.read_csv(native_path)
    stats = annual_metrics(selected, native)
    stats.to_csv(dst/"Fig2_old_vs_new_annual_metrics.csv",index=False,float_format="%.17g")
    merged = selected.merge(native[["date","swe_i_mm","swe_n_mm"]].assign(date=lambda t:pd.to_datetime(t.date)),on="date",suffixes=("_new","_old"))
    merged["change_i_mm"] = merged.swe_i_mm_new-merged.swe_i_mm_old
    merged["change_n_mm"] = merged.swe_n_mm_new-merged.swe_n_mm_old
    merged.to_csv(dst/"Fig2_old_vs_new_daily.csv",index=False,date_format="%Y-%m-%d",float_format="%.17g")
    source_paths = [src/n for n in ("core_watershed_daily.csv","static_noah.npz","core_methods.json","cache_raw_validation.csv")]
    source_paths += [Path(args.extraction_code),Path(args.grid_metadata),native_path]
    evidence = dst/"provenance"
    evidence.mkdir(exist_ok=True)
    for path in source_paths:
        shutil.copy2(path,evidence/path.name)
    receipt = {
        "status":"Reused authoritative archived spatial common-grid extraction; not a fresh raw model extraction",
        "n_daily_records":len(core), "years":list(range(2018,2027)),
        "state_support":"D23 UTC; full WY2018--2025, October--June WY2026",
        "spatial_operator":"iSnobal 50m states area-overlap averaged to native Noah-MP 1km grid by rasterio/GDAL Resampling.average; Noah native cell states unchanged",
        "original_minimum_regrid_valid_coverage":.5,
        "original_coverage_operator":"GDAL average of finite native-state mask over source raster overlap. Source-edge normalization must be audited separately; this is not automatically full-target-cell native coverage.",
        "fixed_footprint":"Intended footprint is every positive fractional watershed-area target cell; archived joint finite coverage is constant 0.9999999403953552 on all dates, within float32 rounding of unity. Spatial caches are required for fresh direct all-cell verification; no daily mask changes are evidenced by the archived coverage.",
        "normalization":"sum(SWE_j * fractional_watershed_area_j)/sum(fractional_watershed_area_j), constant denominator",
        "finite_zero_policy":"Retained as valid snow-free states",
        "daily_common_valid_area_min":float(core.valid_common_area_fraction.min()),
        "daily_common_valid_area_max":float(core.valid_common_area_fraction.max()),
        "watershed_projected_area_m2_float64":float(np.sum(w,dtype=np.float64)*grid_meta["cell_area_m2"]),
        "float32_legacy_area_m2":methods["watershed_area_m2"],
        "coverage_precision_audit":precision,
        "raw_spotcheck_receipt":"18 archived checks, two per WY; iSnobal SWE maxdifference0mm and Noah SWE maxdifference0mm. Spotchecks were not rerun in this revision.",
        "full_spatial_cache_availability":"Original full NPZ model spatial daily caches remain on HPC, absent from local package; no native-grid 1D means were remapped.",
        "target_grid":grid_meta,
        "source_files":[{"path":str(p.resolve()),"sha256":sha256(p)} for p in source_paths],
        "outputs":{p.name:sha256(p) for p in dst.iterdir() if p.is_file()},
    }
    (dst/"model_operator_provenance.json").write_text(json.dumps(receipt,indent=2),encoding="utf-8")
    print(json.dumps({"records":len(core),"footprint_cells":int((w>0).sum()),"common_coverage_min":float(core.valid_common_area_fraction.min()),"output":str(dst)},indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-dir",required=True,help="Authoritative existing actual_figures data directory")
    p.add_argument("--native-daily",required=True,help="Previous native-grid Fig.2 daily CSV, for comparison only")
    p.add_argument("--extraction-code",required=True,help="Archived extract_core.py used to generate source summary")
    p.add_argument("--grid-metadata",required=True,help="Existing grid_diagnostics.json with Noah-MP grid")
    p.add_argument("--output-dir",required=True)
    import_verified(p.parse_args())


if __name__=="__main__":
    main()
