#!/usr/bin/env python3
"""Plot the validated 14-flight ASO SWE elevation-band scores without pooling.

Keep this helper beside Fig5_ASO_Performance.py in code/f5/. It uses that
module's canonical score validation, flight dates, and elevation definitions.
Input is the complete corrected 224-row score CSV; this presentation-only
helper exports exactly its 84 SWE Lower/Middle/Upper rows. It opens no rasters
and runs no models. Public source distributions do not include research data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.patches import Rectangle
from Fig5_ASO_Performance import FLIGHT_DATES, BANDS, load_scores_csv

MODELS = ('iSnobal', 'Noah-MP')
ELEVATION_BANDS = ('Lower', 'Middle', 'Upper')
METRICS = ('mb_mm', 'rmse_mm', 'r')
TITLES = ('MB (mm)', 'RMSE (mm)', 'Spatial r')
BASENAME = 'JHM_S_ASO_SWE_elevation_performance'
WIDTH_IN = 6.5
HEIGHT_IN = 7.2
FONT_PT = 8.05


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def select_band_rows(rows):
    """Validate and return each canonical SWE elevation key exactly once."""
    selected = [r for r in rows if r['variable'] == 'SWE'
                and r['elevation_band'] in ELEVATION_BANDS]
    lookup = {}
    for row in selected:
        key = (row['date'], row['model'], row['elevation_band'])
        if key in lookup:
            raise ValueError(f'Duplicate SWE elevation key: {key}')
        lookup[key] = row
    expected = {(d, m, b) for d in FLIGHT_DATES for m in MODELS
                for b in ELEVATION_BANDS}
    if set(lookup) != expected:
        raise ValueError('Expected exactly all 84 SWE elevation-band rows')
    return [lookup[d, m, b] for d in FLIGHT_DATES for m in MODELS
            for b in ELEVATION_BANDS]


def score_arrays(rows):
    """Map chronological flights by elevation column; undefined r stays NaN."""
    selected = select_band_rows(rows)
    lookup = {(r['date'], r['model'], r['elevation_band']): r
              for r in selected}
    arrays = {}
    for model in MODELS:
        for metric in METRICS:
            values = [[lookup[d, model, b][metric] for b in ELEVATION_BANDS]
                      for d in FLIGHT_DATES]
            arrays[model, metric] = np.array(values, dtype=float)
    return arrays


def metric_scales(arrays):
    """Fixed shared scales refuse out-of-range values instead of clipping."""
    scales = {'mb_mm': (-300., 300.), 'rmse_mm': (0., 400.), 'r': (-1., 1.)}
    for metric, (low, high) in scales.items():
        for model in MODELS:
            vals = arrays[model, metric]
            finite = vals[np.isfinite(vals)]
            if finite.size and (finite.min() < low or finite.max() > high):
                raise ValueError(f'{metric} exceeds shared display scale; update it explicitly')
    return scales


def configure_font(font_path):
    if font_path:
        prop = font_manager.FontProperties(fname=font_path)
        family = prop.get_name()
        # Resolve the requested font file even when another family copy exists.
        font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist
                                           if f.name != family]
        font_manager.fontManager.addfont(str(font_path))
        font_manager.fontManager._findfont_cached.cache_clear()
        actual = Path(font_manager.findfont(font_manager.FontProperties(family=family)))
        # Path.resolve() can retain '..' in the Windows sandbox while the font
        # registry returns the normalized absolute spelling of the same file.
        normalize = lambda p: os.path.normcase(os.path.normpath(os.path.abspath(p)))
        if normalize(actual) != normalize(font_path):
            raise ValueError('Requested font did not resolve to its supplied file')
    else:
        family = 'DejaVu Sans'
    plt.rcParams.update({'font.family': family, 'font.size': FONT_PT,
                         'axes.titlesize': FONT_PT, 'axes.labelsize': FONT_PT,
                         'xtick.labelsize': FONT_PT, 'ytick.labelsize': FONT_PT,
                         'legend.fontsize': FONT_PT, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'axes.linewidth': .5,
                         'savefig.facecolor': 'white'})
    return family


def json_array(array):
    return [[None if not np.isfinite(v) else float(v) for v in row]
            for row in array]


def render_figure(rows, output, font_path=None):
    family = configure_font(font_path)
    arrays = score_arrays(rows)
    scales = metric_scales(arrays)
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(WIDTH_IN, HEIGHT_IN), sharey='row')
    fig.subplots_adjust(left=.163, right=.987, top=.936, bottom=.137,
                        wspace=.21, hspace=.29)
    cmaps = {'mb_mm': 'RdBu_r', 'rmse_mm': 'YlOrRd', 'r': 'BrBG'}
    artists = {}
    panel_audit = []
    for model_index, model in enumerate(MODELS):
        for metric_index, (metric, title) in enumerate(zip(METRICS, TITLES)):
            ax = axes[model_index, metric_index]
            values = arrays[model, metric]
            cmap = plt.get_cmap(cmaps[metric]).copy()
            cmap.set_bad('#dddddd')
            low, high = scales[metric]
            norm = (colors.TwoSlopeNorm(vmin=low, vcenter=0., vmax=high)
                    if metric in ('mb_mm', 'r') else colors.Normalize(low, high))
            image = ax.imshow(np.ma.masked_invalid(values), aspect='auto',
                              interpolation='nearest', cmap=cmap, norm=norm)
            artists[metric] = image
            ax.set_xticks(range(3), ELEVATION_BANDS)
            ax.set_yticks(range(len(FLIGHT_DATES)), FLIGHT_DATES)
            ax.tick_params(axis='both', length=0, pad=3)
            ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
            ax.set_yticks(np.arange(-.5, len(FLIGHT_DATES), 1), minor=True)
            ax.grid(which='minor', color='white', linewidth=.45)
            ax.tick_params(which='minor', length=0)
            for spine in ax.spines.values():
                spine.set_visible(False)
            panel = chr(97 + model_index * 3 + metric_index)
            ax.set_title(f'({panel}) {title}', pad=5, loc='left')
            for row_index, column_index in np.argwhere(~np.isfinite(values)):
                ax.add_patch(Rectangle((column_index - .5, row_index - .5), 1, 1,
                                       facecolor='#dddddd', hatch='///',
                                       edgecolor='#999999', linewidth=0.))
                ax.text(column_index, row_index, 'NA', color='#222222',
                        ha='center', va='center', fontsize=FONT_PT,
                        bbox={'facecolor': '#dddddd', 'edgecolor': 'none', 'pad': .4})
            panel_audit.append({'panel': panel, 'model': model, 'metric': metric,
                                'flight_rows': list(FLIGHT_DATES),
                                'elevation_columns': list(ELEVATION_BANDS),
                                'values': json_array(values),
                                'unavailable_mask': (~np.isfinite(values)).tolist(),
                                'scale': [low, high], 'colormap': cmaps[metric]})
        box = axes[model_index, 0].get_position()
        fig.text(box.x0, box.y1 + .043, model, ha='left', va='bottom',
                 fontsize=FONT_PT)
    for index, metric in enumerate(METRICS):
        box = axes[1, index].get_position()
        cax = fig.add_axes([box.x0, .071, box.width, .016])
        colorbar = fig.colorbar(artists[metric], cax=cax, orientation='horizontal')
        ticks = {'mb_mm': [-300, 0, 300], 'rmse_mm': [0, 200, 400], 'r': [-1, 0, 1]}[metric]
        colorbar.set_ticks(ticks)
        colorbar.ax.tick_params(labelsize=FONT_PT, length=2, width=.5, pad=2)
        colorbar.outline.set_linewidth(.5)
    fig.text(.987, .018, 'Hatched NA: undefined spatial r', ha='right', va='bottom',
             fontsize=FONT_PT)
    fig.canvas.draw()
    # Verify every visible text object at the native final publication size.
    text_objects = [t for t in fig.findobj(matplotlib.text.Text)
                    if t.get_visible() and t.get_text()]
    if min(t.get_fontsize() for t in text_objects) < FONT_PT:
        raise ValueError('A visible label is smaller than the configured print font')
    renderer = fig.canvas.get_renderer()
    clipped = []
    for t in text_objects:
        bbox = t.get_window_extent(renderer)
        bounds = fig.bbox
        if bbox.x0 < bounds.x0 - .1 or bbox.x1 > bounds.x1 + .1 or bbox.y0 < bounds.y0 - .1 or bbox.y1 > bounds.y1 + .1:
            clipped.append(t.get_text())
    if clipped:
        raise ValueError(f'Text falls outside figure canvas: {clipped}')
    for extension in ('pdf', 'png'):
        fig.savefig(output / f'{BASENAME}.{extension}', dpi=600,
                    metadata=({'Creator': 'Supp_ASO_SWE_Elevation.py',
                               'CreationDate': None, 'ModDate': None}
                              if extension == 'pdf' else None))
    fig.savefig(output / 'elevation_actual_size_96dpi.png', dpi=96)
    plt.close(fig)
    return {'width_in': WIDTH_IN, 'height_in': HEIGHT_IN, 'minimum_font_pt': FONT_PT,
            'font_family': family, 'panels': panel_audit,
            'text_objects_checked': len(text_objects), 'clipped_text': clipped}


def export_exact_band_rows(source, destination):
    """Preserve every original scalar string in the 84-row public schema."""
    with source.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        rows = [r for r in reader if r['variable'] == 'SWE'
                and r['elevation_band'] in ELEVATION_BANDS]
    with destination.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scores-csv', required=True, type=Path)
    parser.add_argument('--output-directory', required=True, type=Path)
    parser.add_argument('--font', type=Path)
    args = parser.parse_args()
    scores = load_scores_csv(args.scores_csv)
    selected = select_band_rows(scores)
    out = args.output_directory
    out.mkdir(parents=True, exist_ok=True)
    audit = render_figure(selected, out, args.font)
    exported = out / 'ASO_SWE_elevation_scores.csv'
    if export_exact_band_rows(args.scores_csv, exported) != 84:
        raise ValueError('Unexpected exported row count')
    undefined = [{'date': row['date'], 'model': row['model'],
                  'elevation_band': row['elevation_band'], 'n_cells': row['n_cells']}
                 for row in selected if row['r'] is None]
    manifest = {'mode': 'validated_score_replay_no_recomputation',
                'input_scores': {'file': args.scores_csv.name,
                                 'sha256': sha256(args.scores_csv)},
                'input_rows': len(scores), 'exported_swe_band_rows': len(selected),
                'elevation_bands_m': {b: list(BANDS[b]) for b in ELEVATION_BANDS},
                'row_order': 'chronological fourteen flight dates, earliest at top',
                'model_rows': list(MODELS), 'metric_columns': list(METRICS),
                'operator': 'Unchanged equal-cell, model-specific finite ASO-model pairs; '
                            'positive watershed eligibility; ASO valid coverage >= 0.5; '
                            'zeros and extremes retained; nominal flight-date D23 model states. '
                            'No pooling between dates or bands.',
                'mb_sign': 'model minus ASO',
                'r': 'spatial Pearson correlation on each flight/model/band finite-pair support',
                'na_style': 'gray hatched cell with NA; distinct from a finite zero',
                'undefined_r': undefined,
                'undefined_r_reason_in_verified_inputs': 'constant Noah-MP lower-band model fields',
                'code_sha256': sha256(Path(__file__)),
                'score_validator_sha256': sha256(Path(sys.modules['Fig5_ASO_Performance'].__file__)),
                'font': {'file': args.font.name, 'sha256': sha256(args.font)} if args.font else
                        {'family': 'DejaVu Sans', 'source': 'Matplotlib bundled default'},
                'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                            'matplotlib': matplotlib.__version__},
                'layout_and_plotted_arrays': audit,
                'outputs': [{'file': p.name, 'sha256': sha256(p)} for p in
                            (out / f'{BASENAME}.pdf', out / f'{BASENAME}.png',
                             out / 'elevation_actual_size_96dpi.png', exported)]}
    (out / 'ASO_SWE_elevation_manifest.json').write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
