"""Recompute F6 and timing sensitivities from approved, fixed-footprint daily means.

This program only reads supplied daily means. It never remaps or simulates snow.
The input must contain date, wy, swe_i_mm, swe_n_mm for every October--June day.
Use --upstream-manifest to record the spatial/clock definition used upstream.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

YEARS = tuple(range(2018, 2027))
MODELS = (("i", "iSnobal"), ("n", "Noah-MP"))
FRACTIONS = (0.4, 0.5, 0.6)


def file_record(path):
    path = Path(path)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def validate_declared_footprint(daily):
    """Check explicit upstream declarations; never infer a changing valid-day mask."""
    for field in ("fixed_model_area_fraction", "fixed_model_cells"):
        if field in daily:
            values = daily[field].to_numpy(float)
            if not np.isfinite(values).all() or not (values > 0).all() or not np.equal(values, values[0]).all():
                raise ValueError(f"Invalid or changing fixed footprint declaration: {field}")
    if "state_hour_utc" in daily and not daily.state_hour_utc.eq(23).all():
        raise ValueError("Mixed or unexpected daily state clock")
    if "daily_common_valid_area_fraction" in daily:
        values = daily.daily_common_valid_area_fraction.to_numpy(float)
        fixed = float(daily.fixed_model_area_fraction.iloc[0]) if "fixed_model_area_fraction" in daily else 1.0
        if not np.isfinite(values).all() or (values < fixed-1e-6).any():
            raise ValueError("Incomplete declared common daily footprint")


def validate_block(block, wy):
    block = block.sort_values("date").copy()
    expected = pd.date_range(f"{wy-1}-10-01", f"{wy}-06-30", freq="D")
    if not pd.DatetimeIndex(block.date).equals(expected):
        raise ValueError(f"WY{wy}: duplicate, missing, or unexpected October--June dates")
    values = block[["swe_i_mm", "swe_n_mm"]].to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"WY{wy}: SWE must be finite and nonnegative on every day")
    return block


def state_metrics(dates, values, fraction):
    """First exact maximum, then first strictly later daily value <= own threshold."""
    dates = pd.DatetimeIndex(dates)
    values = np.asarray(values, dtype=float)
    if len(values) != len(dates) or not len(values):
        raise ValueError("Missing or unaligned dates and states")
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("SWE must be finite and nonnegative")
    if not (0 < fraction < 1):
        raise ValueError("Crossing fraction must lie strictly between zero and one")
    p = int(np.argmax(values))
    peak = float(values[p])
    if peak <= 0:
        raise ValueError("Nonpositive peak: own-peak depletion diagnostic is undefined")
    later = np.flatnonzero(values[p+1:] <= fraction * peak)
    c = p+1+int(later[0]) if len(later) else None
    return dict(peak_index=p, peak_swe_mm=peak, peak_date=str(dates[p].date()),
                crossing_threshold_swe_mm=float(fraction*peak),
                crossing_date=str(dates[c].date()) if c is not None else "",
                peak_to_crossing_days=int((dates[c]-dates[p]).days) if c is not None else None,
                crossing_detected_by_jun30=c is not None,
                swe_on_crossing_mm=float(values[c]) if c is not None else None)


def near_metrics(dates, values, fraction):
    values = np.asarray(values, dtype=float)
    dates = pd.DatetimeIndex(dates)
    selected = np.flatnonzero(values >= fraction*values.max())
    components = np.split(selected, np.flatnonzero(np.diff(selected) > 1)+1)
    return dict(first_near_max_date=str(dates[selected[0]].date()),
                last_near_max_date=str(dates[selected[-1]].date()),
                qualifying_days=int(len(selected)),
                first_to_last_span_inclusive_days=int((dates[selected[-1]]-dates[selected[0]]).days)+1,
                contiguous_component_count=len(components),
                longest_component_days=max(map(len, components)))


def calculate(daily):
    if daily.duplicated(["wy", "date"]).any():
        raise ValueError("Duplicate model daily states")
    thresholds, peaks, deltas, near = [], [], [], []
    for wy in YEARS:
        start, end = pd.Timestamp(wy-1, 10, 1), pd.Timestamp(wy, 6, 30)
        block = validate_block(daily.loc[(daily.wy == wy) & daily.date.between(start, end)], wy)
        dates = pd.DatetimeIndex(block.date)
        for fraction in FRACTIONS:
            row = dict(wy=wy, common_start_date=str(start.date()), common_end_date=str(end.date()),
                       daily_pairs=len(block), fraction_of_own_peak=fraction)
            for tag, model in MODELS:
                result = state_metrics(dates, block[f"swe_{tag}_mm"], fraction)
                row.update({f"{tag}_{k}": v for k, v in result.items() if k != "peak_index"})
                if fraction == 0.5:
                    peak_date, half_date = result["peak_date"], result["crossing_date"]
                    peaks.append(dict(wy=wy, model=model, peak_swe_mm=result["peak_swe_mm"],
                                      peak_date=peak_date, half_date=half_date,
                                      peak_day=(pd.Timestamp(peak_date)-start).days+1,
                                      half_day=(pd.Timestamp(half_date)-start).days+1 if half_date else None,
                                      duration_days=result["peak_to_crossing_days"],
                                      valid_days=len(block), required_days=len(block),
                                      status="observed" if half_date else "censored"))
            row["delta_peak_date_n_minus_i_days"] = (pd.Timestamp(row["n_peak_date"])-pd.Timestamp(row["i_peak_date"])).days
            both = row["i_crossing_detected_by_jun30"] and row["n_crossing_detected_by_jun30"]
            row["delta_crossing_date_n_minus_i_days"] = (pd.Timestamp(row["n_crossing_date"])-pd.Timestamp(row["i_crossing_date"])).days if both else None
            row["delta_duration_n_minus_i_days"] = row["n_peak_to_crossing_days"]-row["i_peak_to_crossing_days"] if both else None
            if both and row["delta_crossing_date_n_minus_i_days"] != row["delta_peak_date_n_minus_i_days"]+row["delta_duration_n_minus_i_days"]:
                raise AssertionError("Calendar decomposition failed")
            thresholds.append(row)
            if fraction == 0.5:
                deltas.append(dict(wy=wy, peak_difference_days=row["delta_peak_date_n_minus_i_days"],
                                   duration_difference_days=row["delta_duration_n_minus_i_days"],
                                   half_difference_days=row["delta_crossing_date_n_minus_i_days"],
                                   i_status="observed" if row["i_crossing_detected_by_jun30"] else "censored",
                                   n_status="observed" if row["n_crossing_detected_by_jun30"] else "censored"))
                for anchor in (1.0, 0.99, 0.95):
                    per_model = {tag: near_metrics(dates, block[f"swe_{tag}_mm"], anchor) for tag, _ in MODELS}
                    first_diff = (pd.Timestamp(per_model["n"]["first_near_max_date"])-pd.Timestamp(per_model["i"]["first_near_max_date"])).days
                    last_diff = (pd.Timestamp(per_model["n"]["last_near_max_date"])-pd.Timestamp(per_model["i"]["last_near_max_date"])).days
                    for tag, model in MODELS:
                        nm = per_model[tag]
                        half = pd.Timestamp(row[f"{tag}_crossing_date"]) if row[f"{tag}_crossing_date"] else None
                        nm.update(wy=wy, model=model, near_max_fraction=anchor,
                                  peak_swe_mm=row[f"{tag}_peak_swe_mm"], argmax_date=row[f"{tag}_peak_date"],
                                  half_date=row[f"{tag}_crossing_date"],
                                  first_near_max_to_existing_half_days=(half-pd.Timestamp(nm["first_near_max_date"])).days if half is not None else None,
                                  last_near_max_to_existing_half_days=(half-pd.Timestamp(nm["last_near_max_date"])).days if half is not None else None,
                                  delta_first_near_max_n_minus_i_days=first_diff,
                                  delta_first_anchor_duration_n_minus_i_days=row["delta_crossing_date_n_minus_i_days"]-first_diff if both else None,
                                  delta_last_near_max_n_minus_i_days=last_diff,
                                  delta_last_anchor_duration_n_minus_i_days=row["delta_crossing_date_n_minus_i_days"]-last_diff if both else None,
                                  delta_existing_half_n_minus_i_days=row["delta_crossing_date_n_minus_i_days"])
                        near.append(nm)
    return pd.DataFrame(thresholds), pd.DataFrame(peaks), pd.DataFrame(deltas), pd.DataFrame(near)


def independent_checks(daily, thresholds):
    """Independent scalar loops: no argmax/search helper and no reuse of calculated indices."""
    checks = []
    for r in thresholds.to_dict("records"):
        block = daily.loc[(daily.wy == r["wy"]) & daily.date.between(r["common_start_date"], r["common_end_date"])].sort_values("date")
        for tag, model in MODELS:
            records = [(d, float(v)) for d, v in zip(block.date, block[f"swe_{tag}_mm"])]
            max_value = max(v for _, v in records)
            peak_date = next(d for d, v in records if v == max_value)
            cross = next((d for d, v in records if d > peak_date and v <= r["fraction_of_own_peak"]*max_value), None)
            computed = dict(peak_swe_mm=max_value, peak_date=str(peak_date.date()),
                            crossing_date=str(cross.date()) if cross is not None else "",
                            peak_to_crossing_days=(cross-peak_date).days if cross is not None else None)
            for key, value in computed.items():
                reference = r[f"{tag}_{key}"]
                match = value == reference or (value is None and pd.isna(reference))
                checks.append(dict(wy=r["wy"], model=model, fraction_of_own_peak=r["fraction_of_own_peak"],
                                   field=key, independent=value, output=reference, exact_match=bool(match)))
                if not match:
                    raise AssertionError(f"Independent scalar check differs: {checks[-1]}")
    return pd.DataFrame(checks)


def independent_near_checks(daily, near):
    checks = []
    model_tag = {model: tag for tag, model in MODELS}
    for r in near.to_dict("records"):
        wy, q, tag = r["wy"], r["near_max_fraction"], model_tag[r["model"]]
        b = daily.loc[(daily.wy == wy) & daily.date.between(f"{wy-1}-10-01", f"{wy}-06-30")].sort_values("date")
        records = list(zip(b.date, b[f"swe_{tag}_mm"].astype(float)))
        peak = max(v for _, v in records)
        chosen = [d for d, v in records if v >= q*peak]
        count, longest, run = 0, 0, 0
        previous = None
        for date in chosen:
            if previous is None or (date-previous).days != 1:
                count += 1
                run = 0
            run += 1
            longest = max(longest, run)
            previous = date
        computed = dict(first_near_max_date=str(chosen[0].date()),
                        last_near_max_date=str(chosen[-1].date()),
                        qualifying_days=len(chosen), contiguous_component_count=count,
                        longest_component_days=longest,
                        first_to_last_span_inclusive_days=(chosen[-1]-chosen[0]).days+1)
        for field, value in computed.items():
            match = value == r[field]
            checks.append(dict(wy=wy, model=r["model"], near_max_fraction=q, field=field,
                               independent=value, output=r[field], exact_match=bool(match)))
            if not match: raise AssertionError(f"Independent near-maximum check differs: {checks[-1]}")
    return pd.DataFrame(checks)


def recrossings(daily, thresholds):
    episodes = []
    for r in thresholds.loc[thresholds.wy == 2026].to_dict("records"):
        for tag, model in MODELS:
            if not r[f"{tag}_crossing_date"]:
                continue
            b = daily.loc[(daily.wy == 2026) & daily.date.between(pd.Timestamp(r[f"{tag}_crossing_date"])+pd.Timedelta(days=1), "2026-06-30")].sort_values("date")
            above = b[f"swe_{tag}_mm"].to_numpy() > r[f"{tag}_crossing_threshold_swe_mm"]
            starts = np.flatnonzero(above & ~np.r_[False, above[:-1]])
            ends = np.flatnonzero(above & ~np.r_[above[1:], False])
            for s, e in zip(starts, ends):
                episodes.append(dict(wy=2026, model=model, fraction_of_own_peak=r["fraction_of_own_peak"],
                                     first_crossing_date=r[f"{tag}_crossing_date"],
                                     later_above_threshold_start_date=str(b.date.iloc[s].date()),
                                     later_above_threshold_end_date=str(b.date.iloc[e].date()),
                                     consecutive_above_threshold_days=int(e-s+1)))
    return pd.DataFrame(episodes, columns=["wy", "model", "fraction_of_own_peak", "first_crossing_date", "later_above_threshold_start_date", "later_above_threshold_end_date", "consecutive_above_threshold_days"])


def plot(daily, peaks, deltas, thresholds, output, font):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib import font_manager
    from matplotlib.ticker import MaxNLocator
    font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({"font.family": font_manager.FontProperties(fname=font).get_name(),
                         "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": .6,
                         "lines.linewidth": 1.1, "savefig.facecolor": "white"})
    colors = {"iSnobal": "#0072B2", "Noah-MP": "#D55E00"}
    def clean(ax):
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(width=.5, length=3)
        ax.yaxis.set_major_locator(MaxNLocator(4))
    def letter(ax, k):
        ax.text(.02, .98, f"({chr(97+k)})", transform=ax.transAxes, ha="left", va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": .5})
    def save(fig, name):
        for ext in ("pdf", "png"):
            fig.savefig(output/f"{name}.{ext}", dpi=600, metadata={"Creator": "timing_common_grid.py"} if ext == "pdf" else None)
        plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 5.4))
    fig.subplots_adjust(left=.10, right=.985, bottom=.09, top=.91, wspace=.3, hspace=.45)
    ax = axs[0, 0]; x = np.arange(len(deltas)); ax.axhline(0, color=".5", lw=.6)
    ax.bar(x-.19, deltas.peak_difference_days, .35, color="#8064A2", label="Peak")
    ax.bar(x+.19, deltas.duration_difference_days, .35, color="#65A695", label="Duration")
    ax.plot(x, deltas.half_difference_days, "kD-", ms=3, lw=.8, label="Half date")
    ax.set_xticks(x, [str(v)[2:] for v in deltas.wy]); ax.set_xlabel("Water year (20xx)")
    ax.set_ylabel("Noah-MP − iSnobal (days)"); ax.set_title("Timing components")
    ax.legend(fontsize=7, frameon=False, ncol=1); clean(ax); letter(ax, 0)
    for k, (ax, wy) in enumerate(zip([axs[0,1], axs[1,0], axs[1,1]], [2018, 2022, 2026]), 1):
        rows = peaks.loc[peaks.wy == wy].set_index("model")
        start = pd.to_datetime(rows.peak_date).min()-pd.Timedelta(days=7)
        end = pd.to_datetime(rows.half_date).max()+pd.Timedelta(days=10)
        ax.axhline(.5, color=".5", lw=.6, ls=":")
        f = daily.loc[(daily.wy == wy) & daily.date.between(start, end)].set_index("date")
        for tag, model in MODELS:
            row = rows.loc[model]; peak = float(row.peak_swe_mm); col = f"swe_{tag}_mm"
            ax.plot(f.index, f[col]/peak, color=colors[model], ls="--" if tag == "n" else "-", label=f"{model}: {peak:.1f} mm")
            for dt, marker in ((row.peak_date, "o"), (row.half_date, "D")):
                date = pd.Timestamp(dt)
                ax.plot(date, f.loc[date, col]/peak, marker, color=colors[model], ms=3)
        ax.set_ylim(0, 1.2); ax.set_xlim(start, end)
        ax.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[1, 15]))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.set_title(f"WY{wy}"); ax.set_ylabel("SWE / own peak")
        ax.legend(frameon=False, fontsize=7, loc="upper right"); clean(ax); letter(ax, k)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color=colors[m], ls="--" if m == "Noah-MP" else "-", label=m) for _, m in MODELS], loc="upper center", bbox_to_anchor=(.5, .995), ncol=2, frameon=False)
    save(fig, "JHM_F6_seasonal_timing")
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    fig.subplots_adjust(left=.12, right=.98, bottom=.22, top=.95)
    for fraction, b in thresholds.groupby("fraction_of_own_peak"):
        ax.plot(b.wy, b.delta_crossing_date_n_minus_i_days, "o-", ms=3, label=f"{100*fraction:.0f}%")
    ax.axhline(0, color=".5", lw=.6); ax.set_ylabel("Crossing difference (days)")
    ax.set_xlabel("Water year"); ax.legend(frameon=False, ncol=3); clean(ax)
    save(fig, "JHM_S_timing_sensitivity")


def write_tables(thresholds, near, output):
    def md(date):
        return pd.Timestamp(date).strftime("%b %d") if date else "unresolved"
    def table(name, label, columns, header, rows, caption):
        text = "\\begin{table}[p]\n\\centering\n\\begingroup\n\\setstretch{1.0}\n\\fontsize{10.5}{12.6}\\selectfont\n\\setlength{\\tabcolsep}{4pt}\n\\renewcommand{\\arraystretch}{1.08}\n"
        text += f"\\begin{{tabular}}{{{columns}}}\n\\toprule\n"+" & ".join(header)+" \\\\\n\\midrule\n"
        text += "\n".join(" & ".join(map(str, row))+" \\\\" for row in rows)
        text += f"\n\\bottomrule\n\\end{{tabular}}\n\\endgroup\n\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{table}}\n"
        (output/name).write_text(text, encoding="utf-8")
    common = "daily SWE on the Noah-MP analysis grid, using the same fixed common model-domain footprint and fractional watershed-area weights as Fig.~2"
    rows = []
    for r in thresholds.loc[thresholds.fraction_of_own_peak == .5].to_dict("records"):
        rows.append([r["wy"], md(r["i_peak_date"]), md(r["n_peak_date"]), f'{r["i_peak_swe_mm"]:.1f}', f'{r["n_peak_swe_mm"]:.1f}', f'{r["delta_peak_date_n_minus_i_days"]:+d}'])
    table("JHM_TS_watershed_peak_dates_20261019.tex", "tab:supp_watershed_exact_peaks", "lrrrrr",
          ["WY", r"\shortstack{iSnobal\\peak date}", r"\shortstack{Noah-MP\\peak date}", r"\shortstack{iSnobal\\peak (mm)}", r"\shortstack{Noah-MP\\peak (mm)}", r"$\Delta$ date (d)"], rows,
          f"First exact maxima of {common} within October--June of each water year. Dates fall in the indicated calendar year; "+r"$\Delta$ is Noah-MP minus iSnobal. SWE peaks are shown to 0.1~mm. The same peaks define the three own-peak crossing fractions in Tables~\ref{tab:supp_watershed_crossing_40}--\ref{tab:supp_watershed_crossing_60}. A selected exact daily maximum is a record-based state metric and need not identify a unique physical onset of depletion.")
    for fraction in FRACTIONS:
        percent = int(100*fraction); rows = []
        for r in thresholds.loc[thresholds.fraction_of_own_peak == fraction].to_dict("records"):
            rows.append([r["wy"], md(r["i_crossing_date"]), md(r["n_crossing_date"]), f'{r["i_peak_to_crossing_days"]} / {r["n_peak_to_crossing_days"]}', f'{r["delta_crossing_date_n_minus_i_days"]:+d}'])
        table(f"JHM_TS_timing_{percent}pct_20261019.tex", f"tab:supp_watershed_crossing_{percent}", "lrrrr",
              ["WY", r"\shortstack{iSnobal\\crossing}", r"\shortstack{Noah-MP\\crossing}", r"\shortstack{Peak-to-crossing (d)\\I / N}", r"$\Delta$ date (d)"], rows,
              f"First strictly postpeak daily crossings at or below {percent}"+r"\% of each model's own October--June exact maximum, using "+common+r". Dates occur in the indicated water year; I / N denotes iSnobal / Noah-MP. Durations are calendar days from each model's own peak to its crossing. $\Delta$ date is Noah-MP minus iSnobal, so a negative value means earlier Noah-MP crossing; this is an intermodel difference, not an error against observations. All crossings are detected by 30 June. There is no interpolation or persistence requirement. The exact decomposition $\Delta$ crossing = $\Delta$ peak + $\Delta$ duration holds for each row. Own-peak thresholds represent different absolute water stores and do not measure melt rate.")
    rows = []
    for wy in (2018, 2022, 2026):
        for fraction in (1., .99, .95):
            for _, model in MODELS:
                r = near.loc[(near.wy == wy) & (near.near_max_fraction == fraction) & (near.model == model)].iloc[0]
                rows.append([wy, model, f"{fraction*100:.0f}", md(r.first_near_max_date), md(r.last_near_max_date), r.qualifying_days, r.contiguous_component_count, f"{r.first_near_max_to_existing_half_days} / {r.last_near_max_to_existing_half_days}"])
    table("JHM_TS_near_maximum_anchors_20261019.tex", "tab:supp_near_maximum_anchors", "llrrrrrr",
          ["WY", "Model", r"$q$ (\%)", "First", "Last", "Days", "Episodes", r"\shortstack{First / last\\to half (d)}"], rows,
          "Descriptive near-maximum sensitivity for WY2018, WY2022 and WY2026 using "+common+r". A qualifying daily state is at least $q$ percent of that model's own October--June maximum. First and last bound all qualifying dates, Days counts qualifying states, and Episodes counts separated contiguous groups; a first-to-last span need not be a continuous plateau. The final column gives durations from the first and last qualifying dates to the original 50\%-of-own-exact-maximum crossing, which is held fixed. Unchanged half dates are not an independent test of crossing-fraction robustness. The 99\% and 95\% tolerances are descriptive choices, not confidence limits or new physical peak dates. Anchor-date plus remaining-duration differences reproduce the original half-date difference.")
    # The compact table already integrated into the current main manuscript.
    rows = []
    for wy in (2018, 2022, 2026):
        for fraction in (1., .99, .95):
            r = near.loc[(near.wy == wy) & (near.near_max_fraction == fraction)].iloc[0]
            rows.append([wy, f"{fraction*100:.0f}", r.delta_first_near_max_n_minus_i_days, r.delta_first_anchor_duration_n_minus_i_days, r.delta_last_near_max_n_minus_i_days, r.delta_last_anchor_duration_n_minus_i_days, r.delta_existing_half_n_minus_i_days])
    table("near_max_timing_20261019.tex", "tab:near_max_timing", "lrrrrrr",
          ["WY", r"Peak (\%)", "First date", "Duration", "Last date", "Duration", "Half date"], rows,
          r"Exact and near-maximum timing anchors for the three displayed seasonal cases, using the same fixed common-grid footprint and fractional watershed-area means as Fig.~2. The percentage is the fraction of each model's own October--June peak. Each anchor pair lists Noah-MP-minus-iSnobal anchor-date and remaining-duration differences (days) to the unchanged first 50\% crossing; the two sum to the half-date difference. First and last qualifying dates may span gaps and do not define a continuous plateau.")


class TimingTests(unittest.TestCase):
    dates = pd.date_range("2024-02-27", periods=6)
    def test_first_exact_tie(self):
        r = state_metrics(self.dates, [0, 10, 10, 6, 5, 0], .5)
        self.assertEqual(r["peak_index"], 1)
        self.assertEqual(r["peak_to_crossing_days"], 3)
    def test_zero_and_equal_threshold_retained(self):
        r = state_metrics(self.dates, [10, 5, 0, 2, 1, 0], .5)
        self.assertEqual(r["crossing_date"], "2024-02-28")
    def test_strictly_after_peak_not_before(self):
        r = state_metrics(self.dates, [0, 0, 10, 6, 5, 0], .5)
        self.assertEqual(r["crossing_date"], "2024-03-02")
    def test_censoring(self):
        r = state_metrics(self.dates, [10, 9, 8, 7, 6, 6], .5)
        self.assertFalse(r["crossing_detected_by_jun30"])
        self.assertIsNone(r["peak_to_crossing_days"])
    def test_zero_peak_undefined(self):
        with self.assertRaises(ValueError): state_metrics(self.dates, np.zeros(6), .5)
    def test_missing_swe_fails(self):
        with self.assertRaises(ValueError): state_metrics(self.dates, [10, 9, np.nan, 4, 1, 0], .5)
    def test_negative_swe_fails(self):
        with self.assertRaises(ValueError): state_metrics(self.dates, [10, 9, 8, 4, 1, -1], .5)
    def test_recovery_does_not_move_first_crossing(self):
        r = state_metrics(self.dates, [10, 5, 9, 8, 7, 0], .5)
        self.assertEqual(r["peak_to_crossing_days"], 1)
    def test_near_max_gaps_counted(self):
        r = near_metrics(self.dates, [10, 9.9, 0, 9.9, 0, 9.9], .99)
        self.assertEqual(r["contiguous_component_count"], 3)
        self.assertEqual(r["qualifying_days"], 4)
        self.assertEqual(r["first_to_last_span_inclusive_days"], 6)
    def test_leap_day_calendar(self):
        d = pd.date_range("2023-10-01", "2024-06-30")
        f = pd.DataFrame({"date": d, "swe_i_mm": 1., "swe_n_mm": 1.})
        self.assertEqual(len(validate_block(f, 2024)), 274)
    def test_calendar_gap_fails(self):
        d = pd.date_range("2023-10-01", "2024-06-30").delete(151)
        with self.assertRaises(ValueError): validate_block(pd.DataFrame({"date": d, "swe_i_mm": 1., "swe_n_mm": 1.}), 2024)
    def test_duplicate_day_fails(self):
        d = pd.date_range("2023-10-01", "2024-06-30")
        d = d.insert(5, d[5])
        with self.assertRaises(ValueError): validate_block(pd.DataFrame({"date": d, "swe_i_mm": 1., "swe_n_mm": 1.}), 2024)
    def test_changing_mask_declaration_fails(self):
        with self.assertRaises(ValueError): validate_declared_footprint(pd.DataFrame({"fixed_model_cells": [842, 841]}))
    def test_missing_declared_support_fails(self):
        with self.assertRaises(ValueError): validate_declared_footprint(pd.DataFrame({"fixed_model_area_fraction": [1., 1.], "daily_common_valid_area_fraction": [1., .99]}))
    def test_mixed_clock_fails(self):
        with self.assertRaises(ValueError): validate_declared_footprint(pd.DataFrame({"state_hour_utc": [23, 0]}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--daily", type=Path)
    p.add_argument("--old-daily", type=Path)
    p.add_argument("--old-thresholds", type=Path)
    p.add_argument("--upstream-manifest", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--figures", type=Path)
    p.add_argument("--font", type=Path)
    p.add_argument("--test", action="store_true")
    args = p.parse_args()
    if args.test:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(TimingTests))
        if not result.wasSuccessful(): raise SystemExit(1)
        if not args.daily: return
    if not all((args.daily, args.old_daily, args.output)):
        p.error("--daily --old-daily --output are required for calculation")
    args.output.mkdir(parents=True, exist_ok=True)
    new = pd.read_csv(args.daily, parse_dates=["date"])
    validate_declared_footprint(new)
    old = pd.read_csv(args.old_daily, parse_dates=["date"])
    th, peaks, deltas, near = calculate(new)
    oldth, oldpeaks, olddeltas, oldnear = calculate(old)
    checks = independent_checks(new, th)
    nearchecks = independent_near_checks(new, near)
    frames = {"timing_threshold_metrics.csv": th, "figure06_timing_metrics.csv": peaks,
              "figure06_timing_differences.csv": deltas, "peak_plateau_sensitivity_all_years.csv": near,
              "timing_independent_scalar_checks.csv": checks,
              "near_maximum_independent_scalar_checks.csv": nearchecks,
              "timing_threshold_wy2026_recrossing_episodes.csv": recrossings(new, th),
              "old_native_vs_common_timing_thresholds.csv": oldth.merge(th, on=["wy", "fraction_of_own_peak"], suffixes=("_old_native", "_new_common")),
              "old_native_vs_common_annual_timing.csv": oldpeaks.merge(peaks, on=["wy", "model"], suffixes=("_old_native", "_new_common")),
              "old_native_vs_common_timing_differences.csv": olddeltas.merge(deltas, on="wy", suffixes=("_old_native", "_new_common")),
              "old_native_vs_common_near_maximum.csv": oldnear.merge(near, on=["wy", "model", "near_max_fraction"], suffixes=("_old_native", "_new_common"))}
    parity = []
    if args.old_thresholds:
        golden = pd.read_csv(args.old_thresholds).set_index(["wy", "fraction_of_own_peak"])
        for r in oldth.to_dict("records"):
            for key in ("i_peak_date", "n_peak_date", "i_crossing_date", "n_crossing_date", "i_peak_to_crossing_days", "n_peak_to_crossing_days", "delta_peak_date_n_minus_i_days", "delta_crossing_date_n_minus_i_days", "delta_duration_n_minus_i_days"):
                match = r[key] == golden.loc[(r["wy"], r["fraction_of_own_peak"]), key]
                parity.append(dict(wy=r["wy"], fraction_of_own_peak=r["fraction_of_own_peak"], field=key, exact_match=bool(match)))
                if not match: raise AssertionError(f"Existing native timing parity failed for {r['wy']} {key}")
        frames["old_native_golden_parity_checks.csv"] = pd.DataFrame(parity)
    for name, frame in frames.items():
        frame.to_csv(args.output/name, index=False, float_format="%.17g")
    write_tables(th, near, args.output)
    if args.figures:
        if not args.font: p.error("--font is required when plotting")
        args.figures.mkdir(parents=True, exist_ok=True)
        plot(new, peaks, deltas, th, args.figures, args.font)
    manifest = dict(analysis="Common Noah-MP analysis-grid, fractionally area-weighted watershed timing",
                    common_horizon="October 1 through June 30, complete actual calendars and leap days",
                    peak="First exact maximum in unsmoothed daily curve",
                    depletion="First strictly later daily value <= 40/50/60 percent of its own model-year peak; no interpolation or persistence",
                    censoring="A nonpositive seasonal peak fails; a crossing absent by June30 remains missing rather than being assigned June30",
                    near_maximum="First/last qualifying daily anchors >=100/99/95 percent of own peak, with original first 50-percent crossing held fixed; disjoint episodes counted",
                    sign="Noah-MP minus iSnobal; delta crossing = delta peak + delta duration",
                    inputs=[file_record(x) for x in (args.daily, args.old_daily, args.old_thresholds, args.upstream_manifest, args.font, Path(__file__)) if x],
                    independent_scalar_checks=len(checks), scalar_checks_all_pass=bool(checks.exact_match.all()),
                    independent_near_maximum_scalar_checks=len(nearchecks), near_maximum_checks_all_pass=bool(nearchecks.exact_match.all()),
                    old_native_golden_checks=len(parity), old_native_golden_all_pass=all(x["exact_match"] for x in parity),
                    synthetic_unit_tests_passed=result.testsRun if args.test else "Not run in this invocation",
                    thresholds=len(th), model_year_peaks=len(peaks), near_maximum_model_year_anchor_rows=len(near),
                    all_crossings_detected=bool(th[["i_crossing_detected_by_jun30", "n_crossing_detected_by_jun30"]].to_numpy().all()),
                    outputs=[file_record(x) for name in frames for x in [args.output/name]]+[file_record(x) for x in sorted(args.output.glob("*.tex"))]+([file_record(args.figures/f"{stem}.{ext}") for stem in ["JHM_F6_seasonal_timing", "JHM_S_timing_sensitivity"] for ext in ["pdf", "png"]] if args.figures else []))
    (args.output/"timing_methods_and_validation.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(deltas.to_string(index=False))
    print(th.loc[th.wy == 2026, ["fraction_of_own_peak", "i_peak_date", "n_peak_date", "i_crossing_date", "n_crossing_date", "delta_crossing_date_n_minus_i_days"]].to_string(index=False))
    print(json.dumps({"checks": len(checks), "old_native_parity_checks": len(parity), "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
