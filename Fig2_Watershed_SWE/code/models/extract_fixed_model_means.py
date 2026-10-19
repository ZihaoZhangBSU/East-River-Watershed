"""Compute Fig.2 fixed-footprint means from existing remapped spatial caches.

Example:
python extract_fixed_model_means.py --cache-dir CACHE --grid grid_common.npz \
    --output-dir OUTPUT --native-daily OLD_NATIVE.csv \
    --archived-common-daily CORE_WATERSHED_DAILY.csv

Caches must be the existing production ``isnobal_wyYYYY_noahgrid_daily.npz``
and ``noahmp_wyYYYY_noahgrid_daily.npz`` spatial arrays. iSnobal's native 50m
states were area-overlap averaged to this Noah-MP analysis grid when those
caches were produced. This program never remaps a one-dimensional basin mean.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from aggregate_fig2_common_models import aggregate_fixed, annual_metrics, check_calendar, sha256


def extract(cache_dir, grid_path, output_dir, native_daily=None, archived_common_daily=None):
    cache_dir, grid_path, output_dir = map(Path,(cache_dir,grid_path,output_dir))
    output_dir.mkdir(parents=True,exist_ok=True)
    with np.load(grid_path,allow_pickle=False) as z:
        weights = z["weights"].astype(np.float64)
        fixed = z["fixed_domain_mask"].astype(bool)
        transform = z["transform"].tolist()
        crs = str(z["crs"])
    if not np.array_equal(fixed,weights>0):
        raise ValueError("Expected the documented positive fractional watershed footprint")
    rows, receipts, auxiliary_rows = [], [], []
    for wy in range(2018,2027):
        caches = {}
        for name in ("isnobal","noahmp"):
            path = cache_dir/f"{name}_wy{wy}_noahgrid_daily.npz"
            with np.load(path,allow_pickle=False) as z:
                caches[name] = {k:z[k] for k in z.files}
            receipt = {"path":str(path.resolve()),"sha256":sha256(path),"wy":wy,"model":name,
                       "keys":list(caches[name]),"bytes":path.stat().st_size}
            receipts.append(receipt)
        i,n = caches["isnobal"], caches["noahmp"]
        dates = pd.to_datetime(i["dates"].astype(str))
        if not dates.equals(pd.to_datetime(n["dates"].astype(str))):
            raise ValueError(f"Model date mismatch WY{wy}")
        means_i,means_n = aggregate_fixed(i["swe_mm"],n["swe_mm"],weights,fixed)
        for k,date in enumerate(dates):
            # Independent scalar-cell accumulation verifies vectorized sums.
            for state,mean in ((i["swe_mm"][k],means_i[k]),(n["swe_mm"][k],means_n[k])):
                direct = sum(float(v)*float(w) for v,w in zip(state[fixed],weights[fixed]))/sum(float(v) for v in weights[fixed])
                if not np.isclose(mean,direct,rtol=0,atol=1e-10):
                    raise AssertionError("Independent scalar weighted mean mismatch")
            rows.append({"date":date,"wy":wy,"swe_i_mm":means_i[k],"swe_n_mm":means_n[k],
                         "fixed_model_area_fraction":1.0,"daily_common_valid_area_fraction":1.0,
                         "fixed_model_cells":int(fixed.sum()),"state_hour_utc":23})
        # Ancillary coverage is audited independently; SWE coverage must not
        # be reused as evidence for precipitation or a mass-accounting term.
        for model,cache in caches.items():
            for var in cache:
                if var == "dates":
                    continue
                a=np.asarray(cache[var])
                if a.shape != (len(dates),*weights.shape):
                    raise ValueError(f"Unexpected cache shape {model}:{var}:{a.shape}")
                valid=np.isfinite(a[:,fixed])
                fractions=np.sum(valid*weights[fixed],axis=1)/weights[fixed].sum()
                aux={"wy":wy,"model":model,"variable":var,"n_dates":len(dates),
                     "valid_area_min":float(fractions.min()),"valid_area_max":float(fractions.max()),
                     "nonfinite_retained_cell_days":int((~valid).sum()),
                     "negative_retained_cell_days":int((a[:,fixed]<0).sum()),
                     "finite_zero_retained_cell_days":int((a[:,fixed]==0).sum())}
                auxiliary_rows.append(aux)
        print(f"WY{wy}: {len(dates)} dates; {fixed.sum()} fixed cells; both SWE fields fully finite",flush=True)
    frame=check_calendar(pd.DataFrame(rows))
    frame.to_csv(output_dir/"Fig2_common_grid_daily.csv",index=False,date_format="%Y-%m-%d",float_format="%.17g")
    pd.DataFrame(auxiliary_rows).to_csv(output_dir/"spatial_cache_variable_coverage.csv",index=False,float_format="%.17g")
    archived_comparison=None
    if archived_common_daily:
        old_common=check_calendar(pd.read_csv(archived_common_daily))
        if not old_common.date.equals(frame.date):
            raise ValueError("Archived common-grid date mismatch")
        archived_comparison={tag:float(np.max(np.abs(frame[f"swe_{tag}_mm"].to_numpy()-old_common[f"swe_{tag}_mm"].to_numpy()))) for tag in ("i","n")}
        if max(archived_comparison.values())>1e-8:
            raise ValueError(f"New aggregation differs from authoritative archived common series: {archived_comparison}")
    if native_daily:
        native=pd.read_csv(native_daily)
        annual_metrics(frame,native).to_csv(output_dir/"Fig2_old_vs_new_annual_metrics.csv",index=False,float_format="%.17g")
        native=check_calendar(native)
        daily=frame[["date","wy","swe_i_mm","swe_n_mm"]].merge(native[["date","swe_i_mm","swe_n_mm"]],on="date",suffixes=("_new","_old"))
        daily["change_i_mm"]=daily.swe_i_mm_new-daily.swe_i_mm_old
        daily["change_n_mm"]=daily.swe_n_mm_new-daily.swe_n_mm_old
        daily.to_csv(output_dir/"Fig2_old_vs_new_daily.csv",index=False,date_format="%Y-%m-%d",float_format="%.17g")
    manifest={
        "status":"Fresh aggregation of all existing production remapped spatial daily cache fields; no new model runs",
        "n_daily_records":len(frame),"total_independent_scalar_checks":len(frame)*2,
        "time_support":"Daily state at23:00UTC; fullWY2018--2025, October1--June30WY2026",
        "grid":{"crs":crs,"transform":transform,"shape":list(weights.shape)},
        "fixed_model_cells":int(fixed.sum()),
        "fixed_model_area_m2":float(weights[fixed].sum()*abs(transform[0]*transform[4]-transform[1]*transform[3])),
        "weights":"Fractional watershed intersection area in each Noah-MP target cell. These are NOT the equal-cell weights used for Fig.5 spatial scores.",
        "fixed_mask":"Every positive fractional watershed-area target cell, unchanged for every date and both models; all retained cell-days freshly verified finite and nonnegative. Missing retained SWE aborts rather than renormalizing.",
        "operator":"SWE model cache is existing production area-overlap average of native iSnobal50m SWE to Noah-MP1km cells (GDAL/rasterio Resampling.average, minimumfinitecoverage0.5); nativeNoah cellSWE unchanged; daily fixed fractional area weighted mean applied to both.",
        "finite_zeros":"Included; no snow-present or magnitude filters.",
        "regridding_verification_scope":"This revision verifies all cached target fields and recomputes all basin means. It does not re-read all underlying raw50m model states. Archived18raw-date spotchecks and extraction operator source provide prior raw-to-cache verification.",
        "native_validity_scope":"All retained target-cell SWE is finite. Native finite-subcell coverage was screened at0.5 when caches were built, but per-day native coverage arrays are not saved in these caches; 100percent finite target support is not a claim of100percent native valid-subcell coverage.",
        "source_extent_audit":"See ../aso/isnobal_extent_coverage_audit.json: all842 weighted target cells fully within iSnobal raster extent under the full-target padded coverage operator.",
        "production_cache_directory_input":str(cache_dir.resolve()),
        "cache_builder_source_manifest":"provenance/cache_builder_sources.json",
        "fresh_vs_archived_common_max_abs_mm":archived_comparison,
        "grid_input":{"path":str(grid_path.resolve()),"sha256":sha256(grid_path)},
        "spatial_cache_inputs":receipts,
        "outputs":{p.name:sha256(p) for p in output_dir.iterdir() if p.is_file() and p.name != "model_operator_provenance.json"},
    }
    (output_dir/"model_operator_provenance.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    return frame


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache-dir",required=True)
    p.add_argument("--grid",required=True)
    p.add_argument("--output-dir",required=True)
    p.add_argument("--native-daily")
    p.add_argument("--archived-common-daily")
    a=p.parse_args()
    extract(a.cache_dir,a.grid,a.output_dir,a.native_daily,a.archived_common_daily)


if __name__=="__main__":
    main()
