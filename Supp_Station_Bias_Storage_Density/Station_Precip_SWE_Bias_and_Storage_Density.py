"""Recreate station precipitation/SWE bias and storage/density diagnostics.

Uses frozen daily station/model extracts; performs no model simulations or new
imputation. Research inputs and generated figures are distributed separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

YEARS = tuple(range(2018, 2027))
STATIONS = {380: "Butte", 737: "Schofield Pass", 1141: "Upper Taylor*"}
SOURCES = ("SNOTEL", "iSnobal", "Noah-MP")
COLORS = {"SNOTEL": "#222222", "iSnobal": "#0072B2", "Noah-MP": "#D55E00"}
MARKERS = {"SNOTEL": "o", "iSnobal": "o", "Noah-MP": "s"}
LINESTYLES = {"SNOTEL": "-", "iSnobal": "-", "Noah-MP": "--"}
INPUT_FILES = ("snotel_daily.csv", "isnobal_station_daily.csv",
               "noahmp_station_daily.csv", "snotel_precip_gapfill_log.csv")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def year_dates(wy: int) -> pd.DatetimeIndex:
    end = f"{wy}-06-30" if wy == 2026 else f"{wy}-09-30"
    return pd.date_range(f"{wy-1}-10-01", end, freq="D")


def validate_daily(frame: pd.DataFrame, name: str, fields: list[str]) -> pd.DataFrame:
    required = {"date", "station_id", "water_year", *fields}
    if missing := required - set(frame.columns):
        raise ValueError(f"{name}: missing columns {sorted(missing)}")
    q = frame.copy()
    q["date"] = pd.to_datetime(q.date, errors="raise")
    if q.date.isna().any() or q.date.dt.tz is not None:
        raise ValueError(f"{name}: dates must be finite timezone-naive calendar dates")
    if not q.date.equals(q.date.dt.normalize()):
        raise ValueError(f"{name}: subdaily timestamps are not accepted")
    if set(q.station_id) != set(STATIONS):
        raise ValueError(f"{name}: expected station IDs {list(STATIONS)}")
    if q.duplicated(["station_id", "date"]).any():
        raise ValueError(f"{name}: duplicate station/date")
    wy = q.date.dt.year + (q.date.dt.month >= 10).astype(int)
    if not q.water_year.eq(wy).all():
        raise ValueError(f"{name}: inconsistent water-year labels")
    expected = pd.DatetimeIndex([])
    for year in YEARS:
        expected = expected.append(year_dates(year))
    q = q.sort_values(["station_id", "date"]).reset_index(drop=True)
    for station in STATIONS:
        if not pd.DatetimeIndex(q.loc[q.station_id.eq(station), "date"]).equals(expected):
            raise ValueError(f"{name}: incomplete/extra calendar dates at station {station}")
    values = q[fields].to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"{name}: SWE, depth and precipitation must be finite and nonnegative")
    return q


def load_inputs(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    definitions = {
        "SNOTEL": (INPUT_FILES[0], {"swe_mm": "swe_mm", "snow_depth_mm": "depth_mm",
                                     "precip_gapfilled_mm": "precip_mm"}),
        "iSnobal": (INPUT_FILES[1], {"swe_nearest_mm": "swe_mm",
                                      "snow_depth_nearest_mm": "depth_mm",
                                      "precip_nearest_mm": "precip_mm"}),
        "Noah-MP": (INPUT_FILES[2], {"swe_mm": "swe_mm", "snow_depth_mm": "depth_mm",
                                     "precip_mm": "precip_mm"}),
    }
    raw, standardized = {}, []
    for source, (filename, mapping) in definitions.items():
        q = validate_daily(pd.read_csv(data_dir / filename), source, list(mapping))
        raw[source] = q
        d = q[["date", "station_id", "water_year", *mapping]].rename(columns=mapping)
        d["source"] = source
        standardized.append(d)
    station = raw["SNOTEL"]
    required = {"precip_observed_mm", "precip_gapfill_source", "isnobal_gapfill_mm"}
    if required - set(station):
        raise ValueError("SNOTEL table must retain precipitation gap-fill provenance")
    observed = station.precip_observed_mm.to_numpy(float)
    if np.isinf(observed).any() or (observed[np.isfinite(observed)] < 0).any():
        raise ValueError("Observed precipitation must be finite/nonnegative or explicitly missing")
    filled = station.precip_observed_mm.isna()
    if not station.loc[~filled, "precip_gapfill_source"].eq("observed").all():
        raise ValueError("Unfilled precipitation rows must be labeled observed")
    if not station.loc[filled, "precip_gapfill_source"].eq("iSnobal nearest pixel").all():
        raise ValueError("Missing observed precipitation must identify the archived iSnobal fill")
    if not np.array_equal(station.loc[~filled, "precip_observed_mm"],
                          station.loc[~filled, "precip_gapfilled_mm"]):
        raise ValueError("Gap filling changed an observed precipitation value")
    i_precip = raw["iSnobal"].set_index(["station_id", "date"]).precip_nearest_mm
    filled_rows = station.loc[filled].copy()
    keys = pd.MultiIndex.from_frame(filled_rows[["station_id", "date"]])
    if not np.array_equal(filled_rows.precip_gapfilled_mm.to_numpy(), i_precip.reindex(keys).to_numpy()):
        raise ValueError("Filled precipitation must equal the paired nearest iSnobal value")
    if not np.array_equal(filled_rows.precip_gapfilled_mm, filled_rows.isnobal_gapfill_mm):
        raise ValueError("Archived gap-fill values are inconsistent")
    log = pd.read_csv(data_dir / INPUT_FILES[3], parse_dates=["date"])
    log_fields = {"date", "station_id", "precip_gapfilled_mm", "precip_gapfill_source"}
    if log_fields - set(log) or log.duplicated(["station_id", "date"]).any():
        raise ValueError("Gap-fill log must contain unique dates, IDs, values and source labels")
    actual = filled_rows.set_index(["station_id", "date"]).sort_index()
    logged = log.set_index(["station_id", "date"]).sort_index()
    if not actual.index.equals(logged.index):
        raise ValueError("Gap-fill log dates differ from missing observed precipitation")
    if not np.array_equal(actual.precip_gapfilled_mm, logged.precip_gapfilled_mm):
        raise ValueError("Gap-fill log values differ from the daily table")
    if not actual.precip_gapfill_source.equals(logged.precip_gapfill_source):
        raise ValueError("Gap-fill log source labels differ from the daily table")
    daily = pd.concat(standardized, ignore_index=True)
    fill_index = pd.MultiIndex.from_frame(filled_rows[["station_id", "date"]])
    daily_keys = pd.MultiIndex.from_frame(daily[["station_id", "date"]])
    daily["source_precip_is_gapfilled"] = daily_keys.isin(fill_index) & daily.source.eq("SNOTEL")
    daily["density_accepted"] = (np.isfinite(daily.swe_mm) & np.isfinite(daily.depth_mm)
                                  & daily.depth_mm.ge(100.0))
    daily["density_kg_m3"] = np.nan
    valid = daily.density_accepted
    daily.loc[valid, "density_kg_m3"] = 1000 * daily.loc[valid, "swe_mm"] / daily.loc[valid, "depth_mm"]
    return daily.sort_values(["station_id", "source", "date"]).reset_index(drop=True), log


def paired_bias(reference: np.ndarray, model: np.ndarray) -> dict:
    """All supplied paired days contribute, including zero-snow/zero-precip days."""
    reference, model = np.asarray(reference, float), np.asarray(model, float)
    if reference.ndim != 1 or model.shape != reference.shape or len(reference) == 0:
        raise ValueError("Bias requires equal nonempty one-dimensional paired arrays")
    if not np.isfinite(reference).all() or not np.isfinite(model).all():
        raise ValueError("Bias inputs are incomplete; no days may be silently dropped")
    ref_mean, mod_mean = float(reference.mean()), float(model.mean())
    mb = float((model - reference).mean())
    return {"n_paired_days": len(reference), "reference_mean": ref_mean,
            "model_mean": mod_mean, "reference_sum": float(reference.sum()),
            "model_sum": float(model.sum()), "mean_bias": mb,
            "relative_mean_bias_pct": 100 * mb / ref_mean if ref_mean != 0 else np.nan,
            "relative_bias_status": "defined" if ref_mean != 0 else "zero_reference_mean"}


def own_peak_ratio(frame: pd.DataFrame, wy: int) -> dict:
    """First exact SWE maximum in Oct-Jun; precipitation through that source's peak."""
    q = frame.sort_values("date")
    expected = year_dates(wy)
    if not pd.DatetimeIndex(q.date).equals(expected):
        raise ValueError("Peak ratio requires the complete documented available year")
    if not np.isfinite(q[["swe_mm", "precip_mm"]].to_numpy(float)).all():
        raise ValueError("Peak ratio cannot use missing SWE or precipitation")
    search = q[q.date.le(pd.Timestamp(f"{wy}-06-30"))]
    maximum = float(search.swe_mm.max())
    peak_date = search.loc[search.swe_mm.eq(maximum), "date"].iloc[0]
    accumulation = q[q.date.le(peak_date)]
    if not pd.DatetimeIndex(accumulation.date).equals(pd.date_range(f"{wy-1}-10-01", peak_date)):
        raise ValueError("Precipitation denominator must include every date through own peak")
    precip = float(accumulation.precip_mm.sum())
    return {"peak_search_start": f"{wy-1}-10-01", "peak_search_end": f"{wy}-06-30",
            "peak_date": peak_date.date().isoformat(), "peak_swe_mm": maximum,
            "precip_to_own_peak_mm": precip, "precip_to_peak_n_days": len(accumulation),
            "peak_swe_precip_ratio": maximum / precip if precip != 0 else np.nan,
            "peak_ratio_status": "defined" if precip != 0 else "zero_precipitation",
            "source_precip_gapfill_days_to_peak": int(accumulation.source_precip_is_gapfilled.sum())}


def mean_density(swe_mm: np.ndarray, depth_mm: np.ndarray) -> dict:
    """Mean of daily bulk densities on source-specific finite depth >=100 mm dates."""
    swe, depth = np.asarray(swe_mm, float), np.asarray(depth_mm, float)
    if swe.ndim != 1 or swe.shape != depth.shape:
        raise ValueError("Density requires one-dimensional paired SWE/depth arrays")
    finite = np.isfinite(swe) & np.isfinite(depth)
    if (swe[finite] < 0).any() or (depth[finite] < 0).any():
        raise ValueError("Negative SWE/depth cannot be silently filtered")
    valid = finite & (depth >= 100)
    density = 1000 * swe[valid] / depth[valid]
    return {"density_mean_kg_m3": float(density.mean()) if valid.any() else np.nan,
            "density_n_days": int(valid.sum()), "density_finite_pair_days": int(finite.sum()),
            "density_below_100mm_days": int((finite & (depth < 100)).sum()),
            "density_nonfinite_pair_days": int((~finite).sum()),
            "density_status": "defined" if valid.any() else "no_depth_support"}


def calculate(daily: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    bias, storage, fills = [], [], []
    for station in STATIONS:
        for wy in YEARS:
            group = daily[daily.station_id.eq(station) & daily.water_year.eq(wy)]
            ref = group[group.source.eq("SNOTEL")].sort_values("date")
            common = {"station_id": station, "station_name": STATIONS[station].replace("*", ""),
                      "wy": wy, "window_start": ref.date.min().date().isoformat(),
                      "window_end": ref.date.max().date().isoformat(), "window_n_days": len(ref),
                      "partial_year": wy == 2026}
            filled = ref[ref.source_precip_is_gapfilled]
            fill_record = {**common, "gapfill_n_days": len(filled),
                           "gapfill_sum_mm": float(filled.precip_mm.sum()),
                           "gapfill_dates": ";".join(filled.date.dt.strftime("%Y-%m-%d")),
                           "gapfill_source": "nearest iSnobal precipitation"}
            fills.append(fill_record)
            for source in SOURCES:
                d = group[group.source.eq(source)].sort_values("date")
                if not d.date.reset_index(drop=True).equals(ref.date.reset_index(drop=True)):
                    raise ValueError("Reference/model dates must match exactly")
                storage.append({**common, "source": source, **own_peak_ratio(d, wy),
                                **mean_density(d.swe_mm.to_numpy(), d.depth_mm.to_numpy())})
                if source == "SNOTEL":
                    continue
                for variable, column, units in [("precipitation", "precip_mm", "mm/day"),
                                                 ("SWE", "swe_mm", "mm")]:
                    bias.append({**common, "model": source, "variable": variable, "units": units,
                                 **paired_bias(ref[column].to_numpy(), d[column].to_numpy()),
                                 "reference_precip_gapfill_n_days": len(filled),
                                 "reference_precip_gapfill_sum_mm": float(filled.precip_mm.sum())})
    return pd.DataFrame(bias), pd.DataFrame(storage), pd.DataFrame(fills)


def verify_fig3(bias: pd.DataFrame, path: Path) -> dict:
    scores = pd.read_csv(path)
    reference = scores.set_index(["station_id", "wy", "model"])
    calculated = bias[bias.variable.eq("SWE")].set_index(["station_id", "wy", "model"])
    if len(reference) != 54 or set(reference.index) != set(calculated.index):
        raise ValueError("Fig. 3 score table must contain the same 54 station/year/model cases")
    q = calculated.loc[reference.index]
    differences = abs(q.mean_bias.to_numpy() - reference.mb.to_numpy())
    if not np.allclose(q.mean_bias, reference.mb, rtol=0, atol=1e-10):
        raise AssertionError("SWE biases disagree with the accepted Fig. 3 scores")
    if not np.array_equal(q.n_paired_days, reference.amount_n_days):
        raise AssertionError("SWE bias day counts disagree with Fig. 3")
    for col, target in [("window_start", "amount_window_start"), ("window_end", "amount_window_end")]:
        if not np.array_equal(q[col], reference[target]):
            raise AssertionError("SWE bias period disagrees with Fig. 3")
    return {"file": path.name, "sha256": sha256(path), "checked_cases": 54,
            "maximum_absolute_mb_difference_mm": float(differences.max()), "status": "passed"}


def set_style(font: Path | None) -> str:
    family = "DejaVu Sans"
    if font is not None:
        if not font.is_file():
            raise FileNotFoundError(font)
        font_manager.fontManager.addfont(str(font))
        family = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({"font.family": family, "font.size": 8, "axes.labelsize": 8,
                         "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8, "axes.linewidth": .6, "pdf.fonttype": 42,
                         "ps.fonttype": 42, "lines.linewidth": 1.1, "savefig.facecolor": "white"})
    return family


def prepare_axes(axes: np.ndarray) -> None:
    for index, ax in enumerate(axes.flat):
        ax.axvspan(2025.5, 2026.45, color="#EEEEEE", zorder=0)
        ax.axvline(2025.5, color="0.7", linewidth=.5, linestyle=":", zorder=1)
        ax.set_xlim(2017.6, 2026.45)
        ax.set_xticks(YEARS, [str(wy)[2:] for wy in YEARS])
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=2.5, width=.5)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        ax.text(.015, 1.015, f"({chr(97+index)})", transform=ax.transAxes,
                va="bottom", ha="left", fontsize=8, clip_on=False)
    for ax, name in zip(axes[0], STATIONS.values()):
        ax.set_title(name, pad=8)
    for ax in axes[-1]:
        ax.set_xlabel("Water year")


def plot_series(ax, years, values, source):
    ax.plot(years, values, color=COLORS[source], marker=MARKERS[source],
            linestyle=LINESTYLES[source], markersize=3.3, markeredgewidth=.45,
            markeredgecolor="white", zorder=3)


def legend_handles(sources):
    return [Line2D([], [], color=COLORS[s], marker=MARKERS[s], linestyle=LINESTYLES[s],
                   markersize=3.5, label=s) for s in sources]


def save_pair(fig, output_dir, stem):
    paths = []
    for suffix in ("pdf", "png"):
        path = output_dir / f"{stem}.{suffix}"
        fig.savefig(path, dpi=450, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def render(bias, storage, output_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(4, 3, figsize=(6.5, 7.0), sharex=True, sharey="row")
    fig.subplots_adjust(left=.105, right=.985, top=.90, bottom=.075, hspace=.23, wspace=.16)
    prepare_axes(axes)
    rows = [("precipitation", "mean_bias", "Precipitation MB (mm day$^{-1}$)"),
            ("precipitation", "relative_mean_bias_pct", "Precipitation relative MB (%)"),
            ("SWE", "mean_bias", "SWE MB (mm)"),
            ("SWE", "relative_mean_bias_pct", "SWE relative MB (%)")]
    for row, (variable, field, label) in enumerate(rows):
        for column, station in enumerate(STATIONS):
            ax = axes[row, column]
            ax.axhline(0, color="0.6", linestyle=":", linewidth=.7, zorder=1)
            for source in SOURCES[1:]:
                d = bias[bias.station_id.eq(station) & bias.variable.eq(variable) & bias.model.eq(source)]
                plot_series(ax, d.wy, d[field], source)
        axes[row, 0].set_ylabel(label)
        values = bias.loc[bias.variable.eq(variable), field].to_numpy()
        lower, upper = min(0, float(np.nanmin(values))), max(0, float(np.nanmax(values)))
        margin = max((upper-lower)*.14, .03)
        axes[row, 0].set_ylim(lower-margin, upper+margin)
    fig.legend(handles=legend_handles(SOURCES[1:]), loc="upper center", ncol=2,
               bbox_to_anchor=(.54, .99), frameon=False, columnspacing=2)
    paths = save_pair(fig, output_dir, "JHM_S_station_precip_swe_bias")
    fig, axes = plt.subplots(2, 3, figsize=(6.5, 3.9), sharex=True, sharey="row")
    fig.subplots_adjust(left=.105, right=.985, top=.83, bottom=.12, hspace=.27, wspace=.16)
    prepare_axes(axes)
    for row, (field, label) in enumerate([("peak_swe_precip_ratio", "Peak SWE / precipitation\nto own peak"),
                                         ("density_mean_kg_m3", "Mean density (kg m$^{-3}$)")]):
        for column, station in enumerate(STATIONS):
            ax = axes[row, column]
            for source in SOURCES:
                d = storage[storage.station_id.eq(station) & storage.source.eq(source)]
                plot_series(ax, d.wy, d[field], source)
        axes[row, 0].set_ylabel(label)
        values = storage[field].to_numpy()
        lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
        margin = (hi-lo)*.14
        axes[row, 0].set_ylim(max(0, lo-margin), hi+margin)
    fig.legend(handles=legend_handles(SOURCES), loc="upper center", ncol=3,
               bbox_to_anchor=(.54, .99), frameon=False, columnspacing=1.8)
    return paths + save_pair(fig, output_dir, "JHM_S_station_storage_density")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Directory containing the four frozen baseline CSV files")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--font", type=Path, help="Optional NimbusSans-Regular.ttf for journal typography")
    parser.add_argument("--fig3-score-csv", type=Path,
                        help="Optional accepted full-year Fig. 3 score table for 54-case cross-check")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    daily, gapfill_log = load_inputs(args.data_dir)
    bias, storage, gapfills = calculate(daily)
    crosscheck = verify_fig3(bias, args.fig3_score_csv) if args.fig3_score_csv else {"status": "not_requested"}
    tables = {"station_daily_diagnostic_support.csv": daily,
              "station_annual_precip_swe_bias.csv": bias,
              "station_annual_storage_density.csv": storage,
              "station_precip_gapfill_support.csv": gapfills,
              "station_precip_gapfill_dates.csv": gapfill_log}
    for filename, frame in tables.items():
        frame.to_csv(args.output_dir / filename, index=False, float_format="%.17g")
    family = set_style(args.font)
    plots = render(bias, storage, args.output_dir)
    receipt = {
        "inputs": {name: sha256(args.data_dir / name) for name in INPUT_FILES},
        "code_sha256": sha256(Path(__file__)),
        "font": {"file": args.font.name, "sha256": sha256(args.font)} if args.font else {"family": family},
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
        "row_counts": {"daily": len(daily), "annual_bias": len(bias), "annual_storage_density": len(storage),
                       "annual_gapfill_support": len(gapfills), "gapfill_dates": len(gapfill_log)},
        "fig3_crosscheck": crosscheck,
        "definitions": {
            "annual_support": "All calendar days Oct1-Sep30 WY2018-2025; Oct1-Jun30 WY2026, including zeros/leap days.",
            "bias": "Mean(model-reference) across every paired daily date; precipitation reference is the archived gap-filled SNOTEL series.",
            "relative_bias": "100*mean(model-reference)/mean(reference), not mean daily percentage errors.",
            "peak": "First exact maximum unsmoothed SWE between Oct1 and Jun30, independently for each source.",
            "ratio": "Own peak SWE divided by that source's complete daily precipitation sum Oct1 through its own peak date, inclusive.",
            "density": "Arithmetic mean of 1000*SWE_mm/depth_mm over finite SWE/depth and depth>=100 mm dates, separately by source over the available annual support; no density cap.",
            "gapfill": "Four supplied station precipitation gaps are replaced by paired nearest iSnobal precipitation. No new imputation. Filled reference is not fully independent of iSnobal.",
            "interpretation_limit": "Ratios use different source-specific peak dates; density uses different accepted date sets. Neither proves a snow-retention efficiency or a unique causal process."
        },
        "outputs": {path.name: sha256(path) for path in plots + [args.output_dir / name for name in tables]}
    }
    (args.output_dir / "generation_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "generated", "row_counts": receipt["row_counts"], "fig3_crosscheck": crosscheck}, indent=2))


if __name__ == "__main__":
    main()
