"""Generate Figure 2: fractional watershed-area means on the Noah-MP grid.

Build model input with code/models/extract_fixed_model_means.py and ASO input
with extract_fig2_common_aso.py. Daily models share one fixed target footprint;
each unfiltered ASO flight uses its own retained targets and reported coverage.
Research data and generated outputs are supplied separately from public code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
from matplotlib import font_manager
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

YEARS = tuple(range(2018, 2027))
FLIGHT_DATES = (
    "2018-03-31", "2018-05-24", "2019-04-07", "2019-06-10",
    "2022-04-21", "2022-05-18", "2023-04-01", "2023-05-23",
    "2024-04-03", "2024-05-20", "2025-04-07", "2025-05-20",
    "2026-03-18", "2026-04-05",
)
I_COLOR = "#0072B2"
N_COLOR = "#D55E00"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = pd.read_csv(data_dir / "models/Fig2_common_grid_daily.csv", parse_dates=["date"])
    aso = pd.read_csv(data_dir / "aso/aso_common_grid_flights.csv", parse_dates=["date"])
    model_required = {"date", "wy", "swe_i_mm", "swe_n_mm", "fixed_model_area_fraction",
                      "daily_common_valid_area_fraction", "fixed_model_cells", "state_hour_utc"}
    aso_required = {"date", "wy", "aso_watershed_mean_mm", "aso_valid_area_fraction",
                    "retained_target_watershed_area_m2", "fixed_model_watershed_area_m2"}
    for name, frame, required in [("model", model, model_required), ("ASO", aso, aso_required)]:
        if missing := required - set(frame):
            raise ValueError(f"{name} table lacks required columns: {sorted(missing)}")
        if frame.date.isna().any() or frame.date.duplicated().any():
            raise ValueError(f"{name} dates must be finite and unique")
        if not frame.date.equals(frame.date.dt.normalize()):
            raise ValueError(f"{name} dates must use nominal calendar dates")
        expected_wy = frame.date.dt.year + (frame.date.dt.month >= 10).astype(int)
        if not frame.wy.eq(expected_wy).all():
            raise ValueError(f"{name} water-year labels conflict with dates")
    model = model.sort_values("date").reset_index(drop=True)
    aso = aso.sort_values("date").reset_index(drop=True)
    expected_dates = pd.DatetimeIndex([])
    for wy in YEARS:
        end = f"{wy}-06-30" if wy == 2026 else f"{wy}-09-30"
        expected_dates = expected_dates.append(pd.date_range(f"{wy-1}-10-01", end))
    if not pd.DatetimeIndex(model.date).equals(expected_dates):
        raise ValueError("Model data must cover every displayed day exactly once (including leap days)")
    if tuple(aso.date.dt.strftime("%Y-%m-%d")) != FLIGHT_DATES:
        raise ValueError("ASO data must contain exactly the 14 documented flight dates")
    values = model[["swe_i_mm", "swe_n_mm"]].to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Common-grid model SWE must be finite and nonnegative")
    if not model[["fixed_model_area_fraction", "daily_common_valid_area_fraction"]].eq(1).all().all():
        raise ValueError("Both model curves require complete finite support on the fixed target footprint")
    if not model.fixed_model_cells.eq(842).all() or not model.state_hour_utc.eq(23).all():
        raise ValueError("Expected 842 fixed watershed targets and 2300 UTC model states")
    fields = ["aso_watershed_mean_mm", "aso_valid_area_fraction",
              "retained_target_watershed_area_m2", "fixed_model_watershed_area_m2"]
    if not np.isfinite(aso[fields].to_numpy(float)).all():
        raise ValueError("ASO target means and areas must be finite")
    fraction = aso.aso_valid_area_fraction.to_numpy(float)
    if not ((fraction > 0) & (fraction <= 1 + 1e-10)).all():
        raise ValueError("ASO retained targets must cover >0 and <=100% of the fixed basin area")
    if not (
        (aso.retained_target_watershed_area_m2 > 0) & (aso.fixed_model_watershed_area_m2 > 0)
    ).all():
        raise ValueError("ASO target areas must be positive")
    if not np.allclose(aso.fixed_model_watershed_area_m2, 748679519.6193271,
                       rtol=0, atol=1e-3):
        raise ValueError("ASO denominator must match the fixed Noah-MP watershed area")
    if not np.allclose(fraction, aso.retained_target_watershed_area_m2 /
                       aso.fixed_model_watershed_area_m2, rtol=0, atol=1e-10):
        raise ValueError("ASO coverage must equal retained / fixed target watershed area")
    return model, aso


def model_statistics(model: pd.DataFrame, aso: pd.DataFrame) -> pd.DataFrame:
    records = []
    for wy, q in model.groupby("wy", sort=True):
        records.append({
            "wy": int(wy), "start_date": q.date.min().date().isoformat(),
            "end_date": q.date.max().date().isoformat(), "n_dates": len(q),
            "pearson_r": float(q.swe_i_mm.corr(q.swe_n_mm)),
            "mean_bias_noah_minus_isnobal_mm": float((q.swe_n_mm - q.swe_i_mm).mean()),
            "aso_flight_count": int(aso.wy.eq(wy).sum()),
        })
    return pd.DataFrame(records)


def render(model: pd.DataFrame, aso: pd.DataFrame, output_dir: Path, font: Path | None) -> list[Path]:
    family = "DejaVu Sans"
    if font is not None:
        if not font.is_file():
            raise FileNotFoundError(font)
        font_manager.fontManager.addfont(str(font))
        family = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({
        "font.family": family, "font.size": 8, "axes.labelsize": 8,
        "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": 0.6, "lines.linewidth": 1.1, "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(3, 3, figsize=(6.5, 5.5), sharey=True)
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.075, top=0.91, hspace=0.37, wspace=0.15)
    all_values = np.concatenate([model[["swe_i_mm", "swe_n_mm"]].to_numpy().ravel(),
                                 aso.aso_watershed_mean_mm.to_numpy()])
    lower = min(0.0, math.floor(float(all_values.min()) / 50) * 50)
    upper = max(50.0, math.ceil(float(all_values.max()) / 50) * 50)
    markers = 0
    for i, (ax, wy) in enumerate(zip(axes.flat, YEARS)):
        q, a = model[model.wy.eq(wy)], aso[aso.wy.eq(wy)]
        ax.plot(q.date, q.swe_i_mm, color=I_COLOR)
        ax.plot(q.date, q.swe_n_mm, color=N_COLOR, linestyle="--")
        ax.scatter(a.date, a.aso_watershed_mean_mm, s=13, color="black",
                   edgecolors="white", linewidths=0.3, zorder=5)
        markers += len(a)
        ax.set_xlim(pd.Timestamp(f"{wy-1}-10-01"), pd.Timestamp(f"{wy}-06-30" if wy == 2026 else f"{wy}-09-30"))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        ax.set_title(f"WY{wy}" + (" (Oct–Jun)" if wy == 2026 else ""))
        ax.set_ylim(lower, upper)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=0.5, length=3)
        ax.text(0.02, 0.98, f"({chr(97+i)})", transform=ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.5})
    if markers != 14:
        raise AssertionError("The figure must plot all 14 ASO flights")
    fig.supylabel("Watershed SWE (mm)", x=0.012, fontsize=8)
    fig.legend(handles=[
        Line2D([], [], color=I_COLOR, label="iSnobal"),
        Line2D([], [], color=N_COLOR, linestyle="--", label="Noah-MP"),
        Line2D([], [], color="black", marker="o", linestyle="", markersize=3,
               label="ASO (flight coverage)"),
    ], loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=3, frameon=False)
    paths = []
    for suffix in ("pdf", "png"):
        path = output_dir / f"JHM_F2_watershed_swe.{suffix}"
        fig.savefig(path, dpi=600, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=here / "data")
    ap.add_argument("--output-dir", type=Path, default=here / "figures")
    ap.add_argument("--font", type=Path, default=None,
                    help="Optional local font file; defaults to Matplotlib's bundled sans-serif")
    args = ap.parse_args()
    model, aso = load_inputs(args.data_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = render(model, aso, args.output_dir, args.font)
    stats_path = args.output_dir / "Fig2_model_statistics.csv"
    model_statistics(model, aso).to_csv(stats_path, index=False)
    info = {
        "figure": "Fig. 2 - Noah-MP-grid area-weighted watershed SWE",
        "operation": "Figure replay from verified spatial-cache aggregation and raw unfiltered ASO remapping; no model simulation",
        "python": platform.python_version(), "numpy": np.__version__,
        "pandas": pd.__version__, "matplotlib": matplotlib.__version__,
        "generator_sha256": sha256(Path(__file__)),
        "model_rows": len(model), "aso_flights": len(aso),
        "model_support": "Same fixed 842 positive-area Noah-MP targets on every date, fractional watershed area weights, 2300 UTC states",
        "aso_support": "Area average of valid supplied source pixels (including finite zeros) in each full target; retain target when full-cell source-valid area fraction >=0.5; then fractional watershed area weighting of retained targets",
        "aso_coverage_denominator": "Fixed 842-target fractional watershed area; coverage describes retained target area, not exact native observed polygon area",
        "statistics_window": "Each panel's displayed period; common October-June interannual metrics are in data/models/Fig2_old_vs_new_annual_metrics.csv",
        "comparison_limit": "ASO markers have flight-specific footprints. Direct matched-flight comparisons must use data/aso/matched_flight_basin_means.csv, which samples all three on identical retained targets.",
        "aso_mean_range_mm": [float(aso.aso_watershed_mean_mm.min()), float(aso.aso_watershed_mean_mm.max())],
        "aso_coverage_range": [float(aso.aso_valid_area_fraction.min()), float(aso.aso_valid_area_fraction.max())],
        "inputs_sha256": {p.name: sha256(p) for p in [
            args.data_dir / "models/Fig2_common_grid_daily.csv",
            args.data_dir / "aso/aso_common_grid_flights.csv"] + ([args.font] if args.font else [])},
        "outputs_sha256": {p.name: sha256(p) for p in figures + [stats_path]},
    }
    receipt = args.output_dir / "Fig2_execution.json"
    receipt.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "output_dir": str(args.output_dir),
                      "model_rows": len(model), "aso_flights": len(aso)}, indent=2))


if __name__ == "__main__":
    main()
