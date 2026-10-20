"""Update F8's common-grid curve and S5's timing-selected loss sample.

F8 panel c retains the separately defined native-grid input/retention accounting.
All computation reads previously simulated, supplied common-grid daily caches.
"""
from pathlib import Path
import argparse
import json
import unittest
import numpy as np
import pandas as pd
from timing_common_grid import file_record


def selected_losses(dates, reference, model, mask, start, end):
    dates = pd.DatetimeIndex(dates)
    if not dates.equals(pd.date_range(dates[0], dates[-1])):
        raise ValueError('Loss calculation requires a complete consecutive daily calendar')
    reference, model = np.asarray(reference, dtype=float), np.asarray(model, dtype=float)
    if reference.shape != model.shape or reference.shape[0] != len(dates):
        raise ValueError('SWE arrays/dates are not aligned')
    x, y = reference[:-1]-reference[1:], model[:-1]-model[1:]
    keep = np.isfinite(x) & np.isfinite(y) & (x > 1.) & (y > 1.)
    keep &= np.broadcast_to(np.asarray(mask, dtype=bool), keep.shape)
    keep &= np.asarray((dates[1:] >= pd.Timestamp(start)) & (dates[1:] <= pd.Timestamp(end)))[:, None, None]
    t, r, c = np.where(keep)
    return pd.DataFrame(dict(date=dates[1:][t].strftime('%Y-%m-%d'), row=r, col=c,
                             isnobal_swe_loss_mm_day=x[keep], noahmp_swe_loss_mm_day=y[keep],
                             water_year=2026))


def metrics(points):
    x = points.isnobal_swe_loss_mm_day.to_numpy(float)
    y = points.noahmp_swe_loss_mm_day.to_numpy(float)
    return dict(n_cell_days=len(points), pearson_r=float(np.corrcoef(x, y)[0, 1]),
                mean_noah_minus_isnobal_loss_mm_day=float(np.mean(y-x)),
                mean_isnobal_loss_mm_day=float(np.mean(x)), mean_noah_loss_mm_day=float(np.mean(y)),
                first_selected_date=str(points.date.min()), last_selected_date=str(points.date.max()))


def verify_old_points(recomputed, old):
    checks = []
    if len(recomputed) != len(old):
        raise AssertionError(f'Original sample count mismatch: {len(recomputed)} != {len(old)}')
    for date, group in old.groupby('date'):
        new = recomputed.loc[recomputed.date == date]
        if len(new) != len(group): raise AssertionError(f'Original per-date count differs: {date}')
        columns = ['isnobal_swe_loss_mm_day', 'noahmp_swe_loss_mm_day']
        original = group[columns].to_numpy(float)
        actual = new[columns].to_numpy(float)
        original = original[np.lexsort((original[:, 1], original[:, 0]))]
        actual = actual[np.lexsort((actual[:, 1], actual[:, 0]))]
        err = float(np.max(np.abs(actual-original))) if len(actual) else 0.
        np.testing.assert_allclose(actual, original, rtol=0, atol=1e-7)
        checks.append(dict(date=date, old_cell_days=len(group), recomputed_cell_days=len(new),
                           max_pair_absolute_difference_mm_day=err, paired_multiset_match=True))
    return pd.DataFrame(checks)


def configure(font, f8=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.family': family,
                         'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8.5 if f8 else 8,
                         'xtick.labelsize': 7.5 if f8 else 7, 'ytick.labelsize': 7.5 if f8 else 7,
                         'legend.fontsize': 7.5 if f8 else 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'axes.linewidth': .65 if f8 else .6, 'lines.linewidth': 1.15 if f8 else 1.1,
                         'savefig.facecolor': 'white', 'legend.frameon': False,
                         'mathtext.fontset': 'custom', 'mathtext.rm': family,
                         'mathtext.it': family+':italic', 'mathtext.bf': family,
                         'mathtext.cal': family, 'mathtext.sf': family})
    return plt


def plot_s5(points, out, font):
    plt = configure(font)
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import MaxNLocator
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    fig.subplots_adjust(left=.12, right=.85, top=.97, bottom=.12)
    x, y = points.isnobal_swe_loss_mm_day, points.noahmp_swe_loss_mm_day
    mesh = ax.hexbin(x, y, gridsize=45, mincnt=1, norm=LogNorm(), cmap='viridis', linewidths=0, rasterized=True)
    assert int(mesh.get_array().sum()) == len(points)
    hi = max(x.max(), y.max())*1.02
    ax.plot([0, hi], [0, hi], color='.3', lw=.7)
    ax.set(xlim=(0, hi), ylim=(0, hi), xlabel='iSnobal net SWE loss (mm d$^{-1}$)', ylabel='Noah-MP net SWE loss (mm d$^{-1}$)')
    fig.colorbar(mesh, ax=ax, label='Cell-day count')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(width=.5, length=3); ax.yaxis.set_major_locator(MaxNLocator(4))
    for ext in ['pdf', 'png']: fig.savefig(out/f'JHM_S5_joint_daily_loss.{ext}', dpi=600)
    plt.close(fig)


def plot_f8(daily, footprint, precipitation, account, start_event, end_event, out, font):
    plt = configure(font, f8=True)
    import matplotlib.dates as mdates
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    IC, NC = '#0072B2', '#D55E00'
    def style(ax):
        ax.spines[['top', 'right']].set_visible(False)
        ax.yaxis.grid(True, color='0.85', ls='--', lw=.45, zorder=0); ax.set_axisbelow(True)
    def panel(ax, label, title):
        ax.set_title(f'{label}  {title}', loc='left', pad=7); style(ax)
    start, end, cutoff = pd.Timestamp('2026-02-01'), pd.Timestamp('2026-06-30'), pd.Timestamp('2026-03-15')
    fig = plt.figure(figsize=(6.5, 6.5))
    gs = fig.add_gridspec(3, 1, left=.115, right=.98, bottom=.105, top=.955, hspace=.82, height_ratios=[1, 1, 1.1])
    a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1], sharex=a); c = fig.add_subplot(gs[2])
    daily = daily.loc[(daily.wy == 2026) & daily.date.between(start, end)].sort_values('date')
    a.axvspan(start_event, end_event, color='.91', zorder=0)
    a.plot(daily.date, daily.swe_i_mm, color=IC, label='iSnobal')
    a.plot(daily.date, daily.swe_n_mm, color=NC, label='Noah-MP')
    a.set_ylim(bottom=0); a.set_ylabel('Watershed SWE (mm)')
    a.legend(loc='upper right', ncol=2, handlelength=1.8)
    panel(a, '(a)', 'WY2026 common-grid snow storage')
    footprint = footprint.loc[footprint.date.between(start, end)].sort_values('date')
    categories = ['both_pct', 'i_only_pct', 'n_only_pct', 'neither_pct', 'unknown_pct']
    np.testing.assert_allclose(footprint[categories].sum(axis=1), 100, rtol=0, atol=1e-7)
    assert footprint.unknown_pct.max() < 1e-8
    colors = ['#8064A2', IC, NC, '#D9DDE2']; labels = ['Both', 'iSnobal only', 'Noah-MP only', 'Neither']
    b.stackplot(footprint.date, *[footprint[col].to_numpy(float) for col in categories[:4]], colors=colors, lw=0)
    precipitation = precipitation.loc[precipitation.date.between(start, end)]
    assert (precipitation.precip_n_mm > 1.).all()
    b.scatter(precipitation.date, np.full(len(precipitation), 102.8), marker='v', s=12, color=NC, edgecolor='.15', lw=.35, clip_on=False, zorder=12)
    b.set_ylim(0, 100); b.set_yticks([0, 25, 50, 75, 100]); b.set_ylabel('Watershed footprint (%)')
    panel(b, '(b)', 'Binary snow-state footprint (SWE > 0 mm)')
    handles = [Patch(facecolor=color, label=label) for color, label in zip(colors, labels)]
    handles.append(Line2D([], [], ls='none', marker='v', ms=3.5, color=NC, label='Precipitation > 1 mm day$^{-1}$'))
    b.legend(handles=handles, ncol=3, loc='upper center', bbox_to_anchor=(.5, -.25), borderaxespad=0, columnspacing=1.5, handletextpad=.45, handlelength=1.3)
    for ax in [a, b]:
        ax.set_xlim(start, end); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.axvline(cutoff, color='.2', ls=(0, (4, 3)), lw=.75)
    r = account.loc[account.water_year == 2026].iloc[0]
    si, sn = r.Snowfall_Oct02_Mar15_I_mm, r.Snowfall_Oct02_Mar15_N_mm
    gi, gn = r.SWE_gain_Oct01_Mar15_I_mm, r.SWE_gain_Oct01_Mar15_N_mm
    es, eg, er = sn/si-1, gn/gi-1, (gn/sn)/(gi/si)-1
    np.testing.assert_allclose([es, eg, er], [r.epsilon_accumulation, r.epsilon_total_corrected, r.epsilon_ablation], rtol=0, atol=1e-12)
    assert abs((1+eg)-(1+es)*(1+er)) < 1e-12
    percentages = [100*es, 100*eg, 100*er]
    bars = c.bar(np.arange(3), percentages, .55, color=['#2166AC', '#B2182B', '#B2182B'], edgecolor='white', lw=.4)
    for bar, pct in zip(bars, percentages):
        c.annotate(f'{pct:+.2f}%', (bar.get_x()+bar.get_width()/2, pct), xytext=(0, 4 if pct >= 0 else -4), textcoords='offset points', ha='center', va='bottom' if pct >= 0 else 'top', fontsize=8)
    c.set_xticks(np.arange(3), ['Snowfall input\n197.0 → 205.2 mm', 'Net SWE gain\n146.3 → 117.6 mm', 'Inferred retention G/S\n0.742 → 0.573'])
    c.axhline(0, color='.2', lw=.65); c.set_ylim(-32, 13); c.set_yticks([-30, -20, -10, 0, 10]); c.set_ylabel('Relative difference (%)')
    panel(c, '(c)', 'Native-grid accounting to 15 March')
    for ext in ['pdf', 'png']: fig.savefig(out/f'JHM_F8_WY2026_snow_storage_diagnosis.{ext}', dpi=600)
    plt.close(fig)


class LossTests(unittest.TestCase):
    dates = pd.date_range('2026-02-20', periods=4)
    def test_and_strict_loss_threshold(self):
        i = np.array([10, 8, 7, 5.])[:, None, None]
        n = np.array([10, 9, 7, 5.])[:, None, None]
        q = selected_losses(self.dates, i, n, np.ones((1, 1), bool), '2026-02-20', '2026-02-23')
        self.assertEqual(q.date.tolist(), ['2026-02-23'])
    def test_missing_state_does_not_bridge(self):
        i = np.array([10, np.nan, 5, 2.])[:, None, None]
        n = np.array([10, 8, 5, 2.])[:, None, None]
        q = selected_losses(self.dates, i, n, np.ones((1, 1), bool), '2026-02-20', '2026-02-23')
        self.assertEqual(q.date.tolist(), ['2026-02-23'])
    def test_mask_zero_and_window(self):
        i = np.array([[10, 10], [7, 7], [4, 4], [0, 0.]])[:, None, :]
        q = selected_losses(self.dates, i, i, [[True, False]], '2026-02-22', '2026-02-22')
        self.assertEqual(len(q), 1); self.assertEqual(q.col.tolist(), [0])
    def test_accumulation_not_loss(self):
        i = np.array([0, 3, 6, 9.])[:, None, None]
        self.assertTrue(selected_losses(self.dates, i, i, [[True]], '2026-02-20', '2026-02-23').empty)
    def test_calendar_gap_fails(self):
        with self.assertRaises(ValueError): selected_losses(self.dates.delete(1), np.ones((3,1,1)), np.ones((3,1,1)), [[True]], '2026-02-20', '2026-02-23')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--isnobal-cache', type=Path); p.add_argument('--noah-cache', type=Path)
    p.add_argument('--grid', type=Path); p.add_argument('--daily', type=Path)
    p.add_argument('--timing', type=Path); p.add_argument('--old-points', type=Path)
    p.add_argument('--baseline-data', type=Path); p.add_argument('--font', type=Path)
    p.add_argument('--output', type=Path); p.add_argument('--figures', type=Path)
    p.add_argument('--test', action='store_true')
    a = p.parse_args()
    if a.test:
        tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(LossTests))
        if not tests.wasSuccessful(): raise SystemExit(1)
        if not a.isnobal_cache: return
    if not all([a.isnobal_cache, a.noah_cache, a.grid, a.daily, a.timing, a.old_points, a.baseline_data, a.output]): p.error('Supply all cache, grid, daily, timing, baseline, and output paths')
    with np.load(a.isnobal_cache, allow_pickle=False) as z: idates, iswe = pd.DatetimeIndex(z['dates']), z['swe_mm'].astype(float)
    with np.load(a.noah_cache, allow_pickle=False) as z: ndates, nswe = pd.DatetimeIndex(z['dates']), z['swe_mm'].astype(float)
    assert idates.equals(ndates)
    assert idates.equals(pd.date_range('2025-10-01', '2026-06-30'))
    with np.load(a.grid, allow_pickle=False) as z:
        mask_key = 'fixed_domain_mask' if 'fixed_domain_mask' in z else 'fixed_model_mask'
        mask, weights = z[mask_key].astype(bool), z['weights'].astype(float)
    assert np.array_equal(mask, np.isfinite(weights) & (weights > 0))
    assert np.isfinite(iswe[:, mask]).all() and np.isfinite(nswe[:, mask]).all()
    timing = pd.read_csv(a.timing)
    wy = timing.loc[timing.wy == 2026]
    start, end = pd.to_datetime(wy.peak_date).min(), pd.to_datetime(wy.half_date).max()
    assert start == pd.Timestamp('2026-02-21') and end == pd.Timestamp('2026-05-10')
    old = pd.read_csv(a.old_points)
    before = selected_losses(idates, iswe, nswe, mask, '2026-02-21', '2026-05-09')
    checks = verify_old_points(before, old)
    new = selected_losses(idates, iswe, nswe, mask, start, end)
    added = new.loc[new.date > '2026-05-09'].copy()
    assert len(new) == len(before)+len(added)
    a.output.mkdir(parents=True, exist_ok=True)
    for name, frame in [('S5_joint_daily_loss_points.csv', new), ('S5_added_May10_cell_days.csv', added), ('S5_original_sample_reproduction_checks.csv', checks)]: frame.to_csv(a.output/name, index=False, float_format='%.17g')
    daily = pd.read_csv(a.daily, parse_dates=['date'])
    d = daily.loc[daily.wy == 2026]
    # Same final means, independently formed from the cached cell arrays.
    wi = np.sum(iswe[:,mask]*weights[mask],axis=1)/weights[mask].sum()
    wn = np.sum(nswe[:,mask]*weights[mask],axis=1)/weights[mask].sum()
    np.testing.assert_allclose(wi, d.swe_i_mm, rtol=0, atol=1e-8)
    np.testing.assert_allclose(wn, d.swe_n_mm, rtol=0, atol=1e-8)
    account_path = a.baseline_data/'F8_h1_recomputed.csv'
    footprint_path = a.baseline_data/'F8_wy2026_snow_state_footprints.csv'
    precipitation_path = a.baseline_data/'F8_wy2026_precipitation_events.csv'
    footprint = pd.read_csv(footprint_path, parse_dates=['date'])
    # Audit unchanged panel b against the freshly reloaded spatial caches.
    bchecks = []
    for r in footprint.to_dict('records'):
        k = idates.get_loc(r['date']); ip, npres = iswe[k,mask] > 0, nswe[k,mask] > 0
        sets = {'both_pct': ip & npres, 'i_only_pct': ip & ~npres, 'n_only_pct': ~ip & npres, 'neither_pct': ~ip & ~npres}
        for key, subset in sets.items():
            actual = float(100*weights[mask][subset].sum()/weights[mask].sum())
            error = actual-r[key]
            if abs(error) > 1e-4: raise AssertionError(f'F8 panel b changed: {r["date"]} {key} {error}')
            bchecks.append(dict(date=str(r['date'].date()), field=key, recomputed_percent=actual, original_percent=r[key], absolute_difference_percent=abs(error)))
    pd.DataFrame(bchecks).to_csv(a.output/'F8_panel_b_fresh_cache_checks.csv',index=False,float_format='%.17g')
    if a.figures:
        if not a.font: p.error('--font required for figures')
        a.figures.mkdir(parents=True, exist_ok=True)
        plot_s5(new, a.figures, a.font)
        plot_f8(daily, footprint, pd.read_csv(precipitation_path, parse_dates=['date']), pd.read_csv(account_path), start, end, a.figures, a.font)
    report = dict(status='Final common-grid dependent refresh; no new model simulation',
                  original_window=['2026-02-21','2026-05-09'], revised_window=[str(start.date()),str(end.date())],
                  old_sample_recomputed=metrics(before), original_sample_matches=bool(checks.paired_multiset_match.all()),
                  revised_sample=metrics(new), added_cell_days=len(added), added_dates=sorted(added.date.unique().tolist()),
                  losses='Previous minus current D23 SWE, assigned to ending date; both losses strictly >1 mm/day, finite paired states, positive fractional watershed target cells; equal weight per selected cell-day',
                  F8_panel_a='Replaced with final common-grid, fixed-footprint, fractional-area SWE; shading uses revised common-grid timing window',
                  F8_panel_b='Archived common-grid presence percentages preserved, verified against new spatial cache reload; triangles remain archived >1mm/day Noah precipitation',
                  F8_panel_c='Native-grid original solid input/gain/retention accounting preserved as separately labeled support; gain is not taken from panel a',
                  F8_panel_b_checks=len(bchecks), F8_panel_b_max_difference_percentage_point=max(x['absolute_difference_percent'] for x in bchecks),
                  tests_passed=tests.testsRun if a.test else 'not run in this invocation',
                  inputs=[file_record(x) for x in [a.isnobal_cache,a.noah_cache,a.grid,a.daily,a.timing,a.old_points,account_path,footprint_path,precipitation_path,Path(__file__)]],
                  outputs=[file_record(x) for x in sorted(a.output.glob('S5*.csv'))]+[file_record(a.output/'F8_panel_b_fresh_cache_checks.csv')]+([file_record(x) for x in sorted(a.figures.glob('JHM_F8*'))]+[file_record(x) for x in sorted(a.figures.glob('JHM_S5*'))] if a.figures else []))
    (a.output/'F8_S5_refresh_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k: report[k] for k in ['old_sample_recomputed','revised_sample','added_cell_days','F8_panel_b_max_difference_percentage_point']},indent=2))


if __name__ == '__main__': main()
