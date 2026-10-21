"""Replay the WY2026 station depth supplement in the accepted compact F4 style.

The 5.85 x 3.40-inch figure uses cached daily depth curves and spatial cell
ranges. No extraction, model simulation, gap filling, resampling, timing
calculation, or SWE timing import is performed. Depth units remain millimetres.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
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
BASENAME = "JHM_S_station_depth"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_inputs(daily):
    """Require the accepted calendar and ordered ranges; never repair data."""
    required = {"date", "station_id", "SNOTEL", "iSnobal", "Noah-MP", "lo", "hi"}
    if not required.issubset(daily):
        raise ValueError("Daily input lacks required depth plotting columns")
    if set(daily.station_id.unique()) != set(STATIONS):
        raise ValueError("Depth input must contain exactly the three accepted stations")
    expected = pd.date_range("2025-10-01", "2026-06-30", freq="D")
    for sid in STATIONS:
        frame = daily.loc[daily.station_id.eq(sid)].sort_values("date")
        if not pd.DatetimeIndex(frame.date).equals(expected):
            raise ValueError(f"Station {sid} must have one row per Oct-June date")
        values = frame[["SNOTEL", "iSnobal", "Noah-MP", "lo", "hi"]].to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError(f"Station {sid} contains nonfinite depth or range data")
        if (frame.lo > frame.hi).any() or (frame.iSnobal < frame.lo).any() or (frame.iSnobal > frame.hi).any():
            raise ValueError(f"Station {sid} has inconsistent iSnobal depth range bounds")


def configure_font(font):
    font = Path(font)
    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=font).get_name()
    normalize = lambda p: os.path.normcase(os.path.normpath(os.path.abspath(p)))
    font_manager.fontManager.ttflist = [entry for entry in font_manager.fontManager.ttflist
                                       if entry.name != family or normalize(entry.fname) == normalize(font)]
    font_manager.fontManager._findfont_cached.cache_clear()
    plt.rcParams.update({
        "font.family": family, "font.size": 8, "axes.labelsize": 8,
        "axes.titlesize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": .6, "lines.linewidth": 1.1, "savefig.facecolor": "white",
    })
    if normalize(font_manager.findfont(family)) != normalize(font):
        raise RuntimeError("Requested publication font was not selected")


def build_figure(daily):
    validate_inputs(daily)
    fig, axes = plt.subplots(3, 1, figsize=(PRINT_WIDTH_IN, PRINT_HEIGHT_IN), sharex=True)
    # Geometry, typography and legend follow the accepted compact main Fig. 4.
    fig.subplots_adjust(left=.12, right=.98, top=.89, bottom=.105, hspace=.13)
    for index, (ax, (sid, name)) in enumerate(zip(axes, STATIONS.items())):
        frame = daily.loc[daily.station_id.eq(sid)].sort_values("date")
        ax.axvspan(pd.Timestamp("2026-03-18"), pd.Timestamp("2026-03-22"),
                   color=GOLD, alpha=.30, lw=0)
        ax.fill_between(frame.date, frame.lo, frame.hi, color=COLORS["iSnobal"], alpha=.16, lw=0)
        for model in COLORS:
            ax.plot(frame.date, frame[model], color=COLORS[model],
                    ls="--" if model == "Noah-MP" else "-")
        # The accepted depth supplement has no peak or disappearance markers.
        # Its depth support must not inherit SWE-derived timing definitions.
        ax.set_ylabel("Depth (mm)")
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
    parser.add_argument("--font", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    configure_font(args.font)
    daily = pd.read_csv(args.plot_data, parse_dates=["date"])
    fig = build_figure(daily)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for extension in ["pdf", "png"]:
        destination = args.output_dir / f"{BASENAME}.{extension}"
        fig.savefig(destination, dpi=600, facecolor="white",
                    metadata={"CreationDate": None, "ModDate": None} if extension == "pdf" else None)
        outputs[destination.name] = {"sha256": sha256(destination), "bytes": destination.stat().st_size}
    fig.savefig(args.output_dir / "station_depth_actual_size_96dpi.png", dpi=96, facecolor="white")
    plt.close(fig)
    receipt = {
        "operation": "Cached-depth-table replay; presentation only", "model_rerun": False,
        "raw_extraction": False, "gap_filling": False, "resampling": False,
        "scientific_data_modified": False, "units": "mm",
        "rows": len(daily), "days_per_station": 273,
        "data_start": "2025-10-01", "data_end": "2026-06-30",
        "display_start": "2026-01-01", "display_end": "2026-06-30",
        "print_width_in": PRINT_WIDTH_IN, "print_height_in": PRINT_HEIGHT_IN,
        "font_sizes_pt": 8,
        "placement_contract": "Place at width 5.85in to preserve the exact compact Fig.4 size and 8pt type.",
        "depth_timing_markers": {"count": 0, "reason": "Accepted depth source has none; SWE timing is not transferred."},
        "range": "Daily minimum-maximum of iSnobal depth pixels within the corresponding Noah-MP cell; spatial range, not a confidence interval.",
        "gold_band": {"start_inclusive": "2026-03-18", "end_inclusive": "2026-03-21", "end_exclusive": "2026-03-22"},
        "inputs_sha256": {"plot_data": sha256(args.plot_data), "font": sha256(args.font)},
        "code_sha256": sha256(__file__), "outputs": outputs,
        "versions": {"python": sys.version, "platform": platform.platform(),
                     "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
    }
    (args.output_dir / "station_depth_generation_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
