"""Generate the annual native-grid snowfall-input/retention supporting figure.

Replays the archived accounting table without changing data or running models.
Provide F8_h1_recomputed.csv and the publication font separately.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager


def configure(font):
    family = font_manager.FontProperties(fname=font).get_name()
    font_manager.fontManager.ttflist = [entry for entry in font_manager.fontManager.ttflist
                                      if entry.name != family]
    font_manager.fontManager.addfont(str(font))
    font_manager.fontManager._findfont_cached.cache_clear()
    actual = font_manager.findfont(font_manager.FontProperties(family=family),
                                   fallback_to_default=False)
    if Path(actual).resolve() != Path(font).resolve():
        raise RuntimeError('Publication font fallback detected')
    plt.rcParams.update({
        'font.family': family, 'font.size': 8., 'axes.labelsize': 8.,
        'axes.titlesize': 8.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
        'legend.fontsize': 7.5, 'lines.linewidth': 1.15, 'axes.linewidth': .65,
        'xtick.major.width': .6, 'ytick.major.width': .6,
        'xtick.major.size': 3, 'ytick.major.size': 3,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'mathtext.fontset': 'custom',
        'mathtext.rm': family, 'mathtext.it': family + ':italic',
        'mathtext.bf': family, 'mathtext.cal': family, 'mathtext.sf': family,
        'legend.frameon': False, 'figure.dpi': 120,
    })


def make_figure(account):
    """Preserve the archived bars, line, limits, styles, and nine-year accounting."""
    d = account.sort_values('water_year')
    if d.water_year.tolist() != list(range(2018, 2027)):
        raise ValueError('Expected one accounting row per WY2018-WY2026')
    si = d.Snowfall_Oct02_Mar15_I_mm.to_numpy(float)
    sn = d.Snowfall_Oct02_Mar15_N_mm.to_numpy(float)
    gi = d.SWE_gain_Oct01_Mar15_I_mm.to_numpy(float)
    gn = d.SWE_gain_Oct01_Mar15_N_mm.to_numpy(float)
    es = d.epsilon_accumulation.to_numpy(float)
    er = d.epsilon_ablation.to_numpy(float)
    eg = d.epsilon_total_corrected.to_numpy(float)
    for a, b in [(es, sn / si - 1), (er, (gn / sn) / (gi / si) - 1),
                 (eg, gn / gi - 1), (1 + eg, (1 + es) * (1 + er))]:
        np.testing.assert_allclose(a, b, atol=1e-12, rtol=0)
    fig, ax = plt.subplots(figsize=(6.5, 3.25))
    fig.subplots_adjust(left=.11, right=.985, bottom=.17, top=.97)
    x = np.arange(len(d)); w = .32
    ax.bar(x - w / 2, 100 * es, width=w, color='#8172B3', edgecolor='black',
           lw=.4, label='Snowfall input')
    ax.bar(x + w / 2, 100 * er, width=w, color='#E29234', edgecolor='black',
           lw=.4, hatch='//', label='Inferred retention')
    ax.plot(x, 100 * eg, color='black', marker='o', mfc='white', mew=.8, ms=3.5,
            label='Net SWE gain')
    ax.axhline(0, color='black', lw=.6)
    ax.set_xticks(x, d.water_year.astype(str)); ax.set_ylim(-30, 10)
    ax.set_ylabel('Relative difference (%)'); ax.set_xlabel('Water year')
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, color='0.85', ls='--', lw=.45, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc='lower left', ncol=3, handlelength=1.5, columnspacing=1.1)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--accounting', type=Path, required=True)
    parser.add_argument('--font', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    configure(args.font)
    fig = make_figure(pd.read_csv(args.accounting))
    args.output.mkdir(parents=True, exist_ok=True)
    stem = 'JHM_S_other_year_annual_input_retention'
    try:
        fig.savefig(args.output / (stem + '.pdf'), facecolor='white', metadata={
            'Title': 'S_other_year_annual_input_retention',
            'Creator': 'Portable archived-data publication replay'})
        fig.savefig(args.output / (stem + '.png'), dpi=600, facecolor='white')
    finally:
        plt.close(fig)
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    receipt = {'operation': 'Archived native-grid accounting replay; snowfall input label',
               'matplotlib_version': matplotlib.__version__,
               'accounting_sha256': sha(args.accounting), 'font_sha256': sha(args.font),
               'code_sha256': sha(Path(__file__)),
               'outputs_sha256': {ext: sha(args.output / (stem + '.' + ext)) for ext in ['pdf', 'png']}}
    (args.output / (stem + '_execution.json')).write_text(
        json.dumps(receipt, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
