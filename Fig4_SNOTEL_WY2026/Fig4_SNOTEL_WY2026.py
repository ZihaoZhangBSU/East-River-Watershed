"""Replay Fig. 4 as a short 5.85 x 3.40-inch plot from verified daily tables.

No raw extraction, model rerun, timing recalculation, or scientific resampling.
Place the resulting PDF at 0.90 times the manuscript text width when that width is
6.5 inches; this is 1:1 placement, not another 10 percent reduction.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

STATIONS = {380: "Butte", 737: "Schofield Pass", 1141: "Upper Taylor*"}
COLORS = {"SNOTEL": "#222222", "iSnobal": "#0072B2", "Noah-MP": "#D55E00"}
GOLD = "#e9b94a"
PRINT_WIDTH_IN = 5.85
PRINT_HEIGHT_IN = 3.40


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_inputs(daily, timing):
    """Reject unverified direct inputs without changing accepted data."""
    required = {"date", "station_id", "SNOTEL", "iSnobal", "Noah-MP", "lo", "hi"}
    if not required.issubset(daily):
        raise ValueError("Daily input lacks required plotting columns")
    if set(daily.station_id.unique()) != set(STATIONS):
        raise ValueError("Daily input must contain exactly the three accepted stations")
    expected = pd.date_range("2025-10-01", "2026-06-30", freq="D")
    for sid in STATIONS:
        f = daily.loc[daily.station_id.eq(sid)].sort_values("date")
        if not pd.DatetimeIndex(f.date).equals(expected):
            raise ValueError(f"Station {sid} must have one row per Oct-June date")
        values = f[["SNOTEL", "iSnobal", "Noah-MP", "lo", "hi"]].to_numpy(float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"Station {sid} contains nonfinite or negative SWE")
        if (f.lo > f.hi).any() or (f.iSnobal < f.lo).any() or (f.iSnobal > f.hi).any():
            raise ValueError(f"Station {sid} has inconsistent iSnobal range bounds")
    keys = ["station_id", "wy", "source", "peak_date", "sdd_date"]
    if not set(keys).issubset(timing):
        raise ValueError("Timing input lacks required columns")
    markers = timing.loc[timing.wy.eq(2026) & timing.station_id.isin(STATIONS)
                         & timing.source.isin(COLORS)].copy()
    if len(markers) != 9 or markers.duplicated(["station_id", "source"]).any():
        raise ValueError("Timing input must contain one row per station and source for WY2026")
    for field in ["peak_date", "sdd_date"]:
        dates = pd.to_datetime(markers[field], errors="raise")
        if dates.isna().any() or not dates.between("2026-01-01", "2026-06-30").all():
            raise ValueError("All 18 accepted peak and disappearance markers must be in the display")
    if (pd.to_datetime(markers.sdd_date) <= pd.to_datetime(markers.peak_date)).any():
        raise ValueError("Each snow-disappearance date must follow its source's own peak date")
    return markers


def configure_font(font):
    font = Path(font).resolve()
    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=font).get_name()
    # Multiple installed copies of a family can otherwise select a different font.
    font_manager.fontManager.ttflist = [entry for entry in font_manager.fontManager.ttflist
                                       if entry.name != family or Path(entry.fname).resolve() == font]
    font_manager.fontManager._findfont_cached.cache_clear()
    plt.rcParams.update({
        "font.family": family, "font.size": 8, "axes.labelsize": 8,
        "axes.titlesize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": .6, "lines.linewidth": 1.1, "savefig.facecolor": "white",
    })
    if Path(font_manager.findfont(family)).resolve() != font:
        raise RuntimeError("Requested publication font was not selected")


def build_figure(daily, timing):
    markers = validate_inputs(daily, timing)
    fig, axes = plt.subplots(3, 1, figsize=(PRINT_WIDTH_IN, PRINT_HEIGHT_IN), sharex=True)
    # Short panels retain the 8 pt print type; titles sit inside the snow-free
    # upper-right corner so their rows no longer consume inter-panel space.
    fig.subplots_adjust(left=.12, right=.98, top=.89, bottom=.105, hspace=.13)
    for index, (ax, (sid, name)) in enumerate(zip(axes, STATIONS.items())):
        f = daily.loc[daily.station_id.eq(sid)].sort_values("date")
        ax.axvspan(pd.Timestamp("2026-03-18"), pd.Timestamp("2026-03-22"),
                   color=GOLD, alpha=.30, lw=0)
        ax.fill_between(f.date, f.lo, f.hi, color=COLORS["iSnobal"], alpha=.16, lw=0)
        for model in COLORS:
            ax.plot(f.date, f[model], color=COLORS[model],
                    ls="--" if model == "Noah-MP" else "-")
            t = markers.loc[markers.station_id.eq(sid) & markers.source.eq(model)].iloc[0]
            peak = pd.Timestamp(t.peak_date)
            value = f.loc[f.date.eq(peak), model].iloc[0]
            ax.plot(peak, value, "o", ms=3, color=COLORS[model])
            ax.plot(pd.Timestamp(t.sdd_date), 10, "v", ms=4, color=COLORS[model])
        ax.set_ylabel("SWE (mm)")
        ax.set_title(name, loc="right", y=.91, pad=0, va="top")
        ax.set_ylim(bottom=-8)
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        ax.set_xlim(pd.Timestamp("2026-01-01"), pd.Timestamp("2026-06-30"))
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=.5, length=3)
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.text(.02, .98, f"({chr(97 + index)})", transform=ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": .5})
    handles = [Line2D([], [], color=COLORS[m], ls="--" if m == "Noah-MP" else "-", label=m)
               for m in COLORS]
    handles += [Patch(facecolor=COLORS["iSnobal"], alpha=.16, label="iSnobal cell range"),
                Patch(facecolor=GOLD, alpha=.3, label="18–21 Mar")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .995),
               ncol=5, frameon=False)
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-data", type=Path, required=True)
    parser.add_argument("--timing-data", type=Path, required=True)
    parser.add_argument("--font", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    configure_font(args.font)
    daily = pd.read_csv(args.plot_data, parse_dates=["date"])
    timing = pd.read_csv(args.timing_data)
    fig = build_figure(daily, timing)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for extension in ["pdf", "png"]:
        destination = args.output_dir / f"JHM_F4_snotel_wy2026.{extension}"
        fig.savefig(destination, dpi=600, facecolor="white",
                    metadata={"CreationDate": None, "ModDate": None} if extension == "pdf" else None)
        outputs[destination.name] = {"sha256": sha256(destination), "bytes": destination.stat().st_size}
    plt.close(fig)
    receipt = {
        "operation": "Cached-table figure replay; layout only", "model_rerun": False,
        "raw_extraction": False, "scientific_data_modified": False,
        "print_width_in": PRINT_WIDTH_IN, "print_height_in": PRINT_HEIGHT_IN,
        "placement_contract": "5.85 in = 0.90 textwidth for a 6.5-in text block; place 1:1",
        "font_sizes_pt": {"ticks_and_legend": 8, "axis_labels_panel_letters_station_titles": 8},
        "inputs_sha256": {"plot_data": sha256(args.plot_data), "timing_data": sha256(args.timing_data),
                          "font": sha256(args.font)},
        "code_sha256": sha256(__file__), "outputs": outputs,
        "versions": {"python": sys.version, "platform": platform.platform(),
                     "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
    }
    (args.output_dir / "generation_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
