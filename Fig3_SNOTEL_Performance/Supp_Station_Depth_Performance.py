#!/usr/bin/env python3
"""Generate the station/year snow-depth counterpart of Fig. 3 (three metrics).

Use the unchanged native nearest-cell daily inputs: SNOTEL snow_depth_mm,
iSnobal snow_depth_nearest_mm, and Noah-MP snow_depth_mm, all in mm.
Amount scores use identical complete dates in all three sources: Oct--Sep
WY2018--2025 and actual Oct--Jun WY2026, retaining zeros and leap days.
Peak dates are each source's own first exact depth maximum in Oct--Jun.
Timing error is model date minus SNOTEL date. No disappearance score is used.
Original archived date labels are retained without timezone shifts. This
checks saved daily data; it does not re-extract raw models or sensor records.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
from datetime import date, timedelta
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

STATIONS = {380: 'Butte', 737: 'Schofield Pass', 1141: 'Upper Taylor'}
MODELS = ('iSnobal', 'Noah-MP')
YEARS = tuple(range(2018, 2027))
SOURCE_COLUMNS = {'SNOTEL': 'snow_depth_mm', 'iSnobal': 'snow_depth_nearest_mm',
                  'Noah-MP': 'snow_depth_mm'}
COLORS = {'iSnobal': '#0072B2', 'Noah-MP': '#D55E00'}
PLOT_FIELDS = ('mb_mm', 'rmse_mm', 'peak_error_days')
BASENAME = 'JHM_S_station_depth_performance'
FONT_PT = 8
TICK_PT = 7
LEGEND_PT = 7
FIGSIZE = (6.5, 6.8)
METRIC_FIELDS = ('station_id', 'station_name', 'wy', 'variable', 'model', 'mb_mm', 'rmse_mm',
                 'peak_error_days', 'observed_peak_date', 'model_peak_date',
                 'observed_peak_depth_mm', 'model_peak_depth_mm', 'observed_peak_tie_count',
                 'model_peak_tie_count', 'amount_window_start', 'amount_window_end',
                 'amount_n_days', 'peak_window_start', 'peak_window_end', 'available_record_end')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calendar_days(start, end):
    while start <= end:
        yield start
        start += timedelta(days=1)


def available_end(wy):
    return date(wy, 6, 30) if wy == 2026 else date(wy, 9, 30)


def amount_days(wy):
    return list(calendar_days(date(wy - 1, 10, 1), available_end(wy)))


def load_source(path, source):
    """Require the same complete, finite native daily contract as current F3."""
    column = SOURCE_COLUMNS[source]
    required = {'date', 'water_year', 'station_id', 'station_name', column}
    groups = {}
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not required.issubset(fields):
            raise ValueError(f'{source}: missing/duplicate required headers')
        for line, row in enumerate(reader, 2):
            if None in row or any(row.get(k) is None for k in required):
                raise ValueError(f'{source}: malformed CSV row {line}')
            try:
                station, wy = int(row['station_id']), int(row['water_year'])
                day, value = date.fromisoformat(row['date']), float(row[column])
            except (ValueError, TypeError) as error:
                raise ValueError(f'{source}: invalid date/metadata/depth at row {line}') from error
            if station not in STATIONS or wy not in YEARS or row['station_name'] != STATIONS[station]:
                raise ValueError(f'{source}: invalid station or water year')
            if row['date'] != day.isoformat() or not date(wy - 1, 10, 1) <= day <= available_end(wy):
                raise ValueError(f'{source}: date outside canonical calendar')
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'{source}: depth must be finite and nonnegative')
            daily = groups.setdefault((station, wy), {})
            if day in daily:
                raise ValueError(f'{source}: duplicate station/date')
            daily[day] = value
    if set(groups) != {(s, y) for s in STATIONS for y in YEARS}:
        raise ValueError(f'{source}: expected all 27 station/year groups')
    for (station, wy), daily in groups.items():
        if set(daily) != set(amount_days(wy)):
            raise ValueError(f'{source}: incomplete daily calendar at {station}/WY{wy}')
    return groups


def first_peak(daily, wy):
    """Source-specific first chronological exact Oct--Jun depth maximum."""
    start, end = date(wy - 1, 10, 1), date(wy, 6, 30)
    season = [(day, value) for day, value in sorted(daily.items())
              if start <= day <= end and value is not None and math.isfinite(value)]
    if not season:
        raise ValueError('No finite Oct--Jun depth values')
    maximum = max(value for _, value in season)
    ties = [day for day, value in season if value == maximum]
    return ties[0], maximum, ties


def amount_metrics(observed, modeled):
    """NumPy reduction, independently checked by a separate math.fsum audit."""
    o, m = np.asarray(observed, float), np.asarray(modeled, float)
    if o.ndim != 1 or o.shape != m.shape or o.size == 0 or not np.isfinite(o).all() or not np.isfinite(m).all():
        raise ValueError('Require equal, nonempty, finite common-date series')
    differences = m - o
    return float(np.mean(differences)), float(np.sqrt(np.mean(differences ** 2)))


def recompute(sources):
    """Compute 54 model rows plus 81 independent source peak events."""
    rows, support, events = [], [], []
    for station, name in STATIONS.items():
        for wy in YEARS:
            days = amount_days(wy)
            daily = {source: sources[source][station, wy] for source in SOURCE_COLUMNS}
            for source, values in daily.items():
                if set(values) != set(days) or not all(math.isfinite(values[d]) and values[d] >= 0 for d in days):
                    raise ValueError(f'{source}: require identical complete finite amount dates')
            peaks = {source: first_peak(values, wy) for source, values in daily.items()}
            support.append({'station_id': station, 'station_name': name, 'wy': wy,
                            'common_amount_n_days': len(days), 'amount_window_start': days[0].isoformat(),
                            'amount_window_end': days[-1].isoformat(),
                            'zero_depth_days': {source: sum(values[d] == 0 for d in days)
                                                for source, values in daily.items()}})
            for source, (peak, maximum, ties) in peaks.items():
                events.append({'station_id': station, 'station_name': name, 'wy': wy, 'source': source,
                               'peak_date': peak.isoformat(), 'peak_depth_mm': maximum,
                               'peak_tie_count': len(ties), 'exact_tie_dates': [d.isoformat() for d in ties],
                               'peak_window_start': date(wy - 1, 10, 1).isoformat(),
                               'peak_window_end': date(wy, 6, 30).isoformat()})
                if source == 'SNOTEL':
                    continue
                mb, rmse = amount_metrics([daily['SNOTEL'][d] for d in days],
                                          [daily[source][d] for d in days])
                observed_peak, observed_maximum, observed_ties = peaks['SNOTEL']
                rows.append({'station_id': station, 'station_name': name, 'wy': wy,
                             'variable': 'Depth', 'model': source, 'mb_mm': mb, 'rmse_mm': rmse,
                             'peak_error_days': (peak - observed_peak).days,
                             'observed_peak_date': observed_peak.isoformat(), 'model_peak_date': peak.isoformat(),
                             'observed_peak_depth_mm': observed_maximum, 'model_peak_depth_mm': maximum,
                             'observed_peak_tie_count': len(observed_ties), 'model_peak_tie_count': len(ties),
                             'amount_window_start': days[0].isoformat(), 'amount_window_end': days[-1].isoformat(),
                             'amount_n_days': len(days), 'peak_window_start': date(wy - 1, 10, 1).isoformat(),
                             'peak_window_end': date(wy, 6, 30).isoformat(), 'available_record_end': available_end(wy).isoformat()})
    validate_metrics(rows)
    return rows, {'common_support': support, 'source_peak_events': events}


def validate_metrics(rows):
    keys = [(int(r['station_id']), int(r['wy']), r['model']) for r in rows]
    expected = {(s, y, m) for s in STATIONS for y in YEARS for m in MODELS}
    if len(keys) != 54 or set(keys) != expected:
        raise ValueError('Require exactly 54 unique station/year/model depth rows')
    observed = {}
    for r in rows:
        station, wy = int(r['station_id']), int(r['wy'])
        if r['variable'] != 'Depth' or r['station_name'] != STATIONS[station]:
            raise ValueError('Metric variable/station mismatch')
        start, end, june = date(wy - 1, 10, 1), available_end(wy), date(wy, 6, 30)
        if (r['amount_window_start'], r['amount_window_end'], r['available_record_end']) != (start.isoformat(), end.isoformat(), end.isoformat()):
            raise ValueError('Incorrect historical full-year / WY2026 partial amount window')
        if int(r['amount_n_days']) != len(amount_days(wy)):
            raise ValueError('Incorrect common amount day count')
        if (r['peak_window_start'], r['peak_window_end']) != (start.isoformat(), june.isoformat()):
            raise ValueError('Incorrect Oct--Jun depth peak window')
        for field in ('mb_mm', 'rmse_mm', 'observed_peak_depth_mm', 'model_peak_depth_mm'):
            if not math.isfinite(float(r[field])):
                raise ValueError(f'Nonfinite {field}')
        if float(r['rmse_mm']) < 0 or float(r['rmse_mm']) + 1e-9 < abs(float(r['mb_mm'])):
            raise ValueError('RMSE must be nonnegative and at least absolute MB')
        op, mp = date.fromisoformat(r['observed_peak_date']), date.fromisoformat(r['model_peak_date'])
        error = float(r['peak_error_days'])
        if not start <= op <= june or not start <= mp <= june or error != (mp - op).days:
            raise ValueError('Invalid depth peak date/error')
        if min(int(r['observed_peak_tie_count']), int(r['model_peak_tie_count'])) < 1:
            raise ValueError('Peak exact-tie counts must be positive')
        observed_key = (op, float(r['observed_peak_depth_mm']), int(r['observed_peak_tie_count']))
        if (station, wy) in observed and observed[station, wy] != observed_key:
            raise ValueError('Both models must use the same observed depth peak')
        observed[station, wy] = observed_key
    return rows


def configure_font(font):
    family = 'DejaVu Sans'
    if font is not None:
        # Normalize before FreeType opens the file, since a long relative path
        # containing parent components can exceed Windows' file-name limit.
        font = Path(os.path.normpath(os.path.abspath(font)))
        family = font_manager.FontProperties(fname=font).get_name()
        font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != family]
        font_manager.fontManager.addfont(str(font))
        font_manager.fontManager._findfont_cached.cache_clear()
        actual = font_manager.findfont(font_manager.FontProperties(family=family))
        normalize = lambda p: os.path.normcase(os.path.normpath(os.path.abspath(p)))
        if normalize(actual) != normalize(font):
            raise ValueError('Requested font does not resolve to the supplied file')
    return family


def make_figure(rows, font=None):
    validate_metrics(rows)
    family = configure_font(font)
    style = {'font.family': family, 'font.size': FONT_PT, 'axes.labelsize': FONT_PT,
             'axes.titlesize': FONT_PT, 'xtick.labelsize': TICK_PT, 'ytick.labelsize': TICK_PT,
             'legend.fontsize': LEGEND_PT, 'pdf.fonttype': 42, 'ps.fonttype': 42,
             'axes.linewidth': .6, 'lines.linewidth': 1.1, 'savefig.facecolor': 'white'}
    lookup = {(int(r['station_id']), int(r['wy']), r['model']): r for r in rows}
    audit = []
    with plt.rc_context(style):
        fig, axes = plt.subplots(3, 3, figsize=FIGSIZE, sharey=True)
        fig.subplots_adjust(left=.12, right=.985, bottom=.08, top=.91, wspace=.23, hspace=.24)
        limits = []
        for field in PLOT_FIELDS:
            extent = max(20., math.ceil(max(abs(float(r[field])) for r in rows) / 20.) * 20.)
            limits.append((0., extent * 1.1) if field == 'rmse_mm' else (-extent * 1.1, extent * 1.1))
        for station_row, (station, name) in enumerate(STATIONS.items()):
            for column, field in enumerate(PLOT_FIELDS):
                ax = axes[station_row, column]
                ax.axhspan(7.5, 8.5, color='#e9b94a', alpha=.25, zorder=0)
                ax.axvline(0, color='.6', lw=.6)
                for year_row, wy in enumerate(YEARS):
                    values = [float(lookup[station, wy, model][field]) for model in MODELS]
                    ax.plot(values, [year_row, year_row], color='.65', lw=.7)
                    for value, model, marker in zip(values, MODELS, ('o', 's')):
                        ax.plot(value, year_row, marker, ms=3.3, mec=COLORS[model],
                                mfc=COLORS[model] if model == 'iSnobal' else 'white', mew=.8)
                        audit.append({'panel': chr(97 + station_row * 3 + column), 'station_id': station,
                                      'wy': wy, 'model': model, 'metric': field, 'value': value,
                                      'year_y': year_row, 'marker': marker, 'color': COLORS[model]})
                ax.set_ylim(8.6, -.6)
                ax.set_yticks(range(9), [str(y)[2:] for y in YEARS])
                ax.set_xlim(limits[column])
                ax.xaxis.set_major_locator(MaxNLocator(3))
                ax.spines[['top', 'right']].set_visible(False)
                if station_row == 0:
                    ax.set_title(('Depth bias', 'Depth RMSE', 'Peak error')[column], pad=8)
                if column == 0:
                    ax.set_ylabel(name + ('*' if station == 1141 else '') + '\nWater year (20xx)')
                if station_row == 2:
                    ax.set_xlabel('mm' if column < 2 else 'days')
                ax.text(.02, .98, f'({chr(97 + station_row * 3 + column)})', transform=ax.transAxes,
                        ha='left', va='top', fontsize=FONT_PT,
                        bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .8, 'pad': .5})
        fig.legend(handles=[Line2D([], [], ls='', marker=marker, mec=COLORS[model],
                                  mfc=COLORS[model] if model == 'iSnobal' else 'white', label=model, ms=4)
                            for model, marker in zip(MODELS, ('o', 's'))],
                   loc='upper center', ncol=2, frameon=False)
        fig.canvas.draw()
        # Locators may keep undrawn ticks outside the view interval. Audit only
        # the labels that are actually eligible for drawing on the page.
        text_candidates = list(fig.legends[0].get_texts())
        for ax in axes.flat:
            text_candidates.extend([ax.title, ax.xaxis.label, ax.yaxis.label, *ax.texts])
            for axis in (ax.xaxis, ax.yaxis):
                low, high = sorted(axis.get_view_interval())
                for tick in axis.get_major_ticks():
                    if low <= tick.get_loc() <= high:
                        text_candidates.extend([tick.label1, tick.label2])
        texts = [t for t in text_candidates if t.get_visible() and t.get_text()]
        if min(t.get_fontsize() for t in texts) < min(FONT_PT, TICK_PT, LEGEND_PT):
            raise ValueError('A visible label is below the print font minimum')
        renderer = fig.canvas.get_renderer()
        clipped = [t.get_text() for t in texts if not fig.bbox.contains(*t.get_window_extent(renderer).p0)
                   or not fig.bbox.contains(*t.get_window_extent(renderer).p1)]
        if clipped:
            raise ValueError(f'Text falls outside the figure canvas: {clipped}')
    return fig, {'width_in': FIGSIZE[0], 'height_in': FIGSIZE[1], 'font_family': family,
                 'minimum_font_pt': min(FONT_PT, TICK_PT, LEGEND_PT),
                 'font_hierarchy_pt': {'labels_titles_panel_letters': FONT_PT,
                                       'ticks': TICK_PT, 'legend': LEGEND_PT},
                 'axis_limits_by_metric': dict(zip(PLOT_FIELDS, limits)),
                 'plotted_markers': audit, 'marker_count': len(audit), 'connector_count': 81,
                 'text_objects_checked': len(texts), 'clipped_text': clipped}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snotel', type=Path)
    parser.add_argument('--isnobal', type=Path)
    parser.add_argument('--noahmp', type=Path)
    parser.add_argument('--metrics', type=Path, help='Replay the complete validated 54-row depth metric table')
    parser.add_argument('--output-directory', type=Path, required=True)
    parser.add_argument('--font', type=Path)
    args = parser.parse_args(argv)
    if args.font is not None:
        args.font = Path(os.path.normpath(os.path.abspath(args.font)))
    paths = {'SNOTEL': args.snotel, 'iSnobal': args.isnobal, 'Noah-MP': args.noahmp}
    if args.metrics:
        if any(paths.values()):
            parser.error('--metrics cannot be combined with daily inputs')
        with args.metrics.open(newline='', encoding='utf-8-sig') as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != METRIC_FIELDS:
                raise ValueError('Depth replay table must use the canonical metric schema')
            rows = validate_metrics(list(reader))
        audit = {'mode': 'validated_metric_replay', 'metric_input_sha256': sha256(args.metrics)}
    else:
        if not all(paths.values()):
            parser.error('Provide all three daily inputs, or --metrics')
        sources = {source: load_source(path, source) for source, path in paths.items()}
        rows, audit = recompute(sources)
        audit['mode'] = 'native_daily_metric_calculation'
        audit['inputs'] = [{'source': source, 'file': path.name, 'sha256': sha256(path),
                            'depth_column': SOURCE_COLUMNS[source], 'units': 'mm',
                            'daily_rows': sum(len(d) for d in sources[source].values())}
                           for source, path in paths.items()]
    out = args.output_directory
    destinations = [out / f'{BASENAME}.pdf', out / f'{BASENAME}.png',
                    out / 'station_depth_metrics.csv', out / 'station_depth_execution.json',
                    out / 'station_depth_actual_size_96dpi.png']
    normalize = lambda p: os.path.normcase(os.path.normpath(os.path.abspath(p)))
    input_paths = {normalize(p) for p in [*paths.values(), args.metrics, args.font] if p}
    if any(normalize(p) in input_paths for p in destinations):
        raise ValueError('Outputs must not replace input files')
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'station_depth_metrics.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=METRIC_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    fig, plot_audit = make_figure(rows, args.font)
    try:
        with plt.rc_context({'pdf.fonttype': 42, 'ps.fonttype': 42}):
            for extension in ('pdf', 'png'):
                fig.savefig(out / f'{BASENAME}.{extension}', dpi=600,
                            metadata=({'Creator': 'Supp_Station_Depth_Performance.py',
                                       'CreationDate': None, 'ModDate': None} if extension == 'pdf' else None))
            fig.savefig(out / 'station_depth_actual_size_96dpi.png', dpi=96)
    finally:
        plt.close(fig)
    audit.update({'metric_rows': len(rows), 'metrics': list(PLOT_FIELDS), 'code_sha256': sha256(__file__),
                  'definitions': {'amount': 'Mean(model minus SNOTEL) and RMSE on complete identical '
                                  'finite dates for all three sources; Oct-Sep historically, Oct-Jun WY2026; '
                                  'zeros and leap days retained, no gap filling or padding.',
                                  'peak': 'Each source own first chronological exact depth maximum Oct-Jun.',
                                  'error': 'Model depth peak date minus SNOTEL depth peak date, calendar days.',
                                  'clock': 'Original saved calendar-date labels, without shifts; '
                                  'exact station/model observation clocks are not assumed identical.'},
                  'font': {'file': args.font.name, 'sha256': sha256(args.font)} if args.font else
                          {'family': 'DejaVu Sans', 'source': 'Matplotlib bundled default'},
                  'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                              'matplotlib': matplotlib.__version__}, 'layout_and_artists': plot_audit,
                  'outputs': [{'file': p.name, 'sha256': sha256(p)} for p in destinations
                              if p.is_file() and p.name != 'station_depth_execution.json']})
    (out / 'station_depth_execution.json').write_text(json.dumps(audit, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
