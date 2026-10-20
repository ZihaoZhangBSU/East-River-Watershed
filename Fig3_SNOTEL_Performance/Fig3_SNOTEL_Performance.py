#!/usr/bin/env python3
"""Generate Figure 3: SNOTEL performance of iSnobal and Noah-MP, WY2018-2026.

Provide the separately supplied 54-row metrics CSV, or build it first with
recompute_fig3_metrics.py. No research data or machine-specific paths are bundled.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

STATIONS = {380: "Butte", 737: "Schofield Pass", 1141: "Upper Taylor*"}
MODELS = ("iSnobal", "Noah-MP")
YEARS = tuple(range(2018, 2027))
COLORS = {"iSnobal": "#0072B2", "Noah-MP": "#D55E00"}
FIELDS = ("mb", "rmse", "peak_error_days", "sdd_error_days")
REQUIRED = {
    "station_id", "station_name", "wy", "model", "mb", "rmse",
    "peak_date", "sdd_date", "peak_error_days", "sdd_error_days",
    "sdd_censored", "available_record_end",
}


def _boolean(value, field):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str) and value in ("True", "False", "true", "false"):
        return value.lower() == "true"
    raise ValueError(f"{field}: expected True or False")


def _date(value, field, optional=False):
    if pd.isna(value) or value == "":
        if optional:
            return pd.NaT
        raise ValueError(f"{field}: missing date")
    try:
        parsed = pd.Timestamp(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"{field}: invalid date") from error
    if pd.isna(parsed) or parsed.tzinfo is not None or parsed != parsed.normalize():
        raise ValueError(f"{field}: expected a calendar date")
    if isinstance(value, str) and value != parsed.strftime("%Y-%m-%d"):
        raise ValueError(f"{field}: expected YYYY-MM-DD")
    return parsed


def validate_metrics(table: pd.DataFrame) -> pd.DataFrame:
    """Require the complete figure grid and internally consistent metric metadata.

    This is validation of a saved metric table, not raw station/model validation.
    Original tables lacking observed_sdd_censored imply an uncensored SNOTEL
    event. New recomputations record that flag explicitly.
    """
    if not table.columns.is_unique or not REQUIRED.issubset(table.columns):
        raise ValueError("Missing or duplicate required metrics columns")
    q = table.copy()
    for field in ("station_id", "wy", *FIELDS):
        try:
            q[field] = pd.to_numeric(q[field].replace("", np.nan), errors="raise")
        except (ValueError, TypeError) as error:
            raise ValueError(f"{field}: invalid number") from error
        values = q[field].to_numpy(dtype=float)
        if field == "sdd_error_days":
            if np.isinf(values).any():
                raise ValueError("SDD errors cannot be infinite")
        elif not np.isfinite(values).all():
            raise ValueError(f"{field}: values must be finite")
        finite = values[np.isfinite(values)]
        if field in ("station_id", "wy", "peak_error_days", "sdd_error_days"):
            if not np.equal(finite, np.floor(finite)).all():
                raise ValueError(f"{field}: values must be whole numbers")
    for field in ("station_id", "wy"):
        q[field] = q[field].astype(int)
    if (q.peak_error_days.abs() > 273).any() or (q.sdd_error_days.abs() > 365).any():
        raise ValueError("Timing error exceeds the supported calendar span")
    keys = list(q[["station_id", "wy", "model"]].itertuples(index=False, name=None))
    expected = set(itertools.product(STATIONS, YEARS, MODELS))
    if len(keys) != len(expected) or set(keys) != expected:
        raise ValueError("Expected exactly 54 unique station/year/model rows")
    names = q.station_id.map({sid: name.rstrip("*") for sid, name in STATIONS.items()})
    if not q.station_name.eq(names).all():
        raise ValueError("Station ID and name do not agree")
    if (q.rmse < 0).any() or (q.rmse + 1e-9 < q.mb.abs()).any():
        raise ValueError("RMSE must be nonnegative and at least the absolute bias")
    if "observed_sdd_censored" not in q:
        q["observed_sdd_censored"] = False
    for field in ("sdd_censored", "observed_sdd_censored"):
        q[field] = q[field].map(lambda value: _boolean(value, field))
    for field in ("peak_date", "sdd_date", "available_record_end"):
        q[field] = q[field].map(lambda value: _date(value, field, field == "sdd_date"))
    for row in q.itertuples(index=False):
        start = pd.Timestamp(row.wy - 1, 10, 1)
        june = pd.Timestamp(row.wy, 6, 30)
        end = june if row.wy == 2026 else pd.Timestamp(row.wy, 9, 30)
        if row.available_record_end != end:
            raise ValueError("Incorrect available record endpoint")
        observed_peak = row.peak_date - pd.Timedelta(days=int(row.peak_error_days))
        if not (start <= row.peak_date <= june and start <= observed_peak <= june):
            raise ValueError("Model and implied observed peaks must be in October-June")
        if row.sdd_censored != pd.isna(row.sdd_date):
            raise ValueError("Model SDD date and censor flag disagree")
        if not row.sdd_censored and not row.peak_date < row.sdd_date <= end - pd.Timedelta(days=6):
            raise ValueError("SDD must follow its own peak and allow seven-day confirmation")
        missing_error = pd.isna(row.sdd_error_days)
        if missing_error != (row.sdd_censored or row.observed_sdd_censored):
            raise ValueError("SDD error must be missing exactly when either event is censored")
        if not missing_error:
            observed_sdd = row.sdd_date - pd.Timedelta(days=int(row.sdd_error_days))
            if not observed_peak < observed_sdd <= end - pd.Timedelta(days=6):
                raise ValueError("Implied observed SDD is outside its confirmation interval")
    for _, group in q.groupby(["station_id", "wy"]):
        if group.observed_sdd_censored.nunique() != 1:
            raise ValueError("Observed censor flag must agree between model rows")
        observed_peaks = group.peak_date - pd.to_timedelta(group.peak_error_days, unit="D")
        if observed_peaks.nunique() != 1:
            raise ValueError("Implied observed peak must agree between model rows")
        finite = group[group.sdd_error_days.notna()]
        observed_sdds = finite.sdd_date - pd.to_timedelta(finite.sdd_error_days, unit="D")
        if observed_sdds.nunique() > 1:
            raise ValueError("Implied observed SDD must agree between model rows")
    return q.sort_values(["station_id", "wy", "model"]).reset_index(drop=True)


def load_metrics(path: Path) -> pd.DataFrame:
    """Read an explicit input path, rejecting duplicate CSV header names."""
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        headers = next(reader, [])
        if len(headers) != len(set(headers)):
            raise ValueError("Duplicate CSV headers")
        if any(len(row) != len(headers) for row in reader):
            raise ValueError("Malformed CSV row length")
    return validate_metrics(pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig"))


def make_figure(table: pd.DataFrame, font: Path | None = None):
    """Plot the reviewed 12-panel figure, preserving its labels and all values."""
    q = validate_metrics(table)
    family = "DejaVu Sans"
    if font is not None:
        font_manager.fontManager.addfont(str(font))
        family = font_manager.FontProperties(fname=font).get_name()
    style = {
        "font.family": family, "font.size": 8, "axes.labelsize": 8,
        "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": .6, "lines.linewidth": 1.1, "savefig.facecolor": "white",
    }
    with plt.rc_context(style):
        fig, axs = plt.subplots(3, 4, figsize=(6.5, 6.8), sharey=True)
        fig.subplots_adjust(left=.12, right=.985, bottom=.08, top=.91, wspace=.23, hspace=.24)
        titles = ("SWE bias", "SWE RMSE", "Peak error", "SDD error")
        limits = []
        for field in FIELDS:
            finite = q.loc[np.isfinite(q[field]), field].abs()
            extent = np.ceil(finite.max() / 20) * 20 if len(finite) else 20
            extent = max(20, extent)  # Also support a valid all-zero/censored synthetic table.
            limits.append((0, extent * 1.1) if field == "rmse" else (-extent * 1.1, extent * 1.1))
        for r, (station_id, name) in enumerate(STATIONS.items()):
            for c, field in enumerate(FIELDS):
                ax = axs[r, c]
                ax.axhspan(7.5, 8.5, color="#e9b94a", alpha=.25, zorder=0)
                ax.axvline(0, color=".6", lw=.6)
                for j, wy in enumerate(YEARS):
                    rows = q[(q.station_id == station_id) & (q.wy == wy)].set_index("model")
                    values = [rows.loc[model, field] for model in MODELS]
                    ax.plot(values, [j, j], color=".65", lw=.7)
                    for value, model, marker in zip(values, MODELS, ("o", "s")):
                        if np.isfinite(value):
                            ax.plot(value, j, marker, ms=3.3, mec=COLORS[model],
                                    mfc=COLORS[model] if model == "iSnobal" else "white", mew=.8)
                        else:
                            ax.text(.98, j, "C", transform=ax.get_yaxis_transform(), ha="right", fontsize=7)
                ax.set_ylim(8.6, -.6)
                ax.set_yticks(range(9), [str(year)[2:] for year in YEARS])
                ax.set_xlim(limits[c])
                ax.xaxis.set_major_locator(MaxNLocator(3))
                ax.spines[["top", "right"]].set_visible(False)
                if r == 0:
                    ax.set_title(titles[c], pad=8)
                if c == 0:
                    ax.set_ylabel(name + "\nWater year (20xx)")
                if r == 2:
                    ax.set_xlabel("mm" if c < 2 else "days")
                ax.text(.02, .98, f"({chr(97 + r * 4 + c)})", transform=ax.transAxes,
                        ha="left", va="top", fontsize=8,
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": .5})
        fig.legend(handles=[
            Line2D([], [], ls="", marker=marker, mec=COLORS[model],
                   mfc=COLORS[model] if model == "iSnobal" else "white", label=model, ms=4)
            for model, marker in zip(MODELS, ("o", "s"))
        ], loc="upper center", ncol=2, frameon=False)
    return fig


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True, help="54-row SNOTEL_figure_data.csv")
    parser.add_argument("--output-dir", type=Path, default=Path("figures"))
    parser.add_argument("--font", type=Path, help="Optional local font for publication appearance")
    args = parser.parse_args(argv)
    table = load_metrics(args.metrics)
    fig = make_figure(table, args.font)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    try:
        for extension in ("pdf", "png"):
            target = args.output_dir / f"F3_snotel_performance.{extension}"
            with plt.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
                fig.savefig(target, dpi=600, facecolor="white")
            outputs[target.name] = sha256(target)
    finally:
        plt.close(fig)
    audit = {
        "figure": "Figure 3: SNOTEL performance", "metric_rows": len(table),
        "plotted_metric_count": len(table) * len(FIELDS),
        "censored_sdd_errors": int(table.sdd_error_days.isna().sum()),
        "metrics_sha256": sha256(args.metrics), "code_sha256": sha256(__file__),
        "font_sha256": sha256(args.font) if args.font else None,
        "font": "supplied font" if args.font else "Matplotlib bundled DejaVu Sans",
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
        "outputs_sha256": outputs,
    }
    (args.output_dir / "fig3_execution.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(f"Figure 3 written to {args.output_dir}; validated {len(table)} metric rows.")


if __name__ == "__main__":
    main()
