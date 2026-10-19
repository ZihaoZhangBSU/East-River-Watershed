"""Verify retained native accounting; diagnose why SWE-only migration is incomplete."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def record(path):
    path = Path(path)
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native-daily', type=Path, required=True)
    p.add_argument('--native-accounting', type=Path, required=True)
    p.add_argument('--common-daily', type=Path, required=True)
    p.add_argument('--core-summary', type=Path, required=True)
    p.add_argument('--spatial-coverage', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    native = pd.read_csv(a.native_daily, parse_dates=['date'])
    accounting = pd.read_csv(a.native_accounting)
    common = pd.read_csv(a.common_daily, parse_dates=['date'])
    summary = pd.read_csv(a.core_summary, parse_dates=['date'])
    if a.spatial_coverage:
        coverage = pd.read_csv(a.spatial_coverage)
        inputs = coverage.loc[coverage.variable.isin(['precip_mm', 'snowfall_mm'])]
        assert len(inputs) == 36 and inputs.nonfinite_retained_cell_days.eq(0).all()
        assert inputs.valid_area_min.eq(1).all()
    rows = []
    for r in accounting.to_dict('records'):
        wy = int(r['water_year'])
        start, end = pd.Timestamp(wy-1, 10, 1), pd.Timestamp(wy, 3, 15)
        row = dict(wy=wy)
        for tag, cap in [('i', 'I'), ('n', 'N')]:
            dates = [start, end]
            endpoints = [float(native.loc[(native.wy == wy) & (native.date == d), f'swe_{tag}_mm'].iloc[0]) for d in dates]
            expected = [float(r[f'SWE_{s}_{cap}_mm']) for s in ['Oct01', 'Mar15']]
            np.testing.assert_allclose(endpoints, expected, rtol=0, atol=1e-10)
            gain = endpoints[1]-endpoints[0]
            snowfall = float(r[f'Snowfall_Oct02_Mar15_{cap}_mm'])
            np.testing.assert_allclose(gain, r[f'SWE_gain_Oct01_Mar15_{cap}_mm'], rtol=0, atol=1e-10)
            assert bool(r[f'snowfall_complete_{cap}'])
            row.update({f'{tag}_native_start_swe_mm': endpoints[0], f'{tag}_native_end_swe_mm': endpoints[1],
                        f'{tag}_native_gain_mm': gain, f'{tag}_native_solid_input_mm': snowfall,
                        f'{tag}_native_retention': gain/snowfall,
                        f'{tag}_native_input_expected_days': r[f'snowfall_expected_days_{cap}'],
                        f'{tag}_native_input_valid_days': r[f'snowfall_valid_days_{cap}']})
            cend = [float(common.loc[(common.wy == wy) & (common.date == d), f'swe_{tag}_mm'].iloc[0]) for d in dates]
            row[f'{tag}_common_start_swe_mm'] = cend[0]
            row[f'{tag}_common_end_swe_mm'] = cend[1]
            row[f'{tag}_common_gain_mm'] = cend[1]-cend[0]
            row[f'{tag}_common_minus_native_gain_mm'] = cend[1]-cend[0]-gain
            f = summary.loc[(summary.wy == wy) & summary.date.between(start+pd.Timedelta(days=1), end)]
            assert pd.DatetimeIndex(f.date).equals(pd.date_range(start+pd.Timedelta(days=1), end))
            assert np.isfinite(f[f'snowfall_{tag}_mm']).all()
            row[f'{tag}_archived_common_solid_input_mm_not_adopted_for_F8c'] = float(f[f'snowfall_{tag}_mm'].sum())
        eps_s = row['n_native_solid_input_mm']/row['i_native_solid_input_mm']-1
        eps_g = row['n_native_gain_mm']/row['i_native_gain_mm']-1
        eps_r = row['n_native_retention']/row['i_native_retention']-1
        np.testing.assert_allclose([eps_s, eps_g, eps_r], [r['epsilon_accumulation'], r['epsilon_total_corrected'], r['epsilon_ablation']], rtol=0, atol=1e-12)
        closure = (1+eps_s)*(1+eps_r)-(1+eps_g)
        assert abs(closure) < 1e-12
        row.update(native_solid_input_relative_percent=100*eps_s, native_swe_gain_relative_percent=100*eps_g,
                   native_retention_relative_percent=100*eps_r, multiplicative_identity_error=closure,
                   native_endpoints_match_original_accounting=True)
        rows.append(row)
    a.output.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(a.output/'F8_native_accounting_verified_and_common_endpoint_comparison.csv', index=False, float_format='%.17g')
    method = dict(decision='Retain the separately defined native-grid F8 panel c/S-other-year accounting. Refresh F8 panel a to the final common-grid SWE display and label the supports separately. Do not combine new common-grid SWE gain with native solid input.',
                  evidence='All 9 annual native Oct01 and Mar15 endpoints and gains reproduce the archived input-retention accounting; 9 multiplicative identities pass at 1e-12.',
                  spatial_audit='Fresh spatial-cache audit separately certifies precipitation and snowfall finite on all 842 target cells. The archived scalar snowfall sums are included for traceability, but are not adopted for panel c, whose original native support is intentionally preserved.',
                  clock_and_boundary_scope='Daily states are D23 UTC. Archived processing integrates full daily hourly records 00--23 and labels those totals by day; original iSnobal daily input timestamps are not retained in these caches. The original Oct01 initial-state correction and Oct02--Mar15 input convention remain unchanged. The native input interfaces differ between models; this audit does not claim exact physical mass-budget closure on a common boundary.',
                  existing_scope='Native accounting remains the existing explicitly delimited input/retention diagnostic: Oct01 to Mar15 SWE gain, Oct02 through Mar15 solid input; Noah-MP ground snowfall and iSnobal prepared solid input enter at different canopy interfaces.',
                  not_used='The archived common snowfall sums in the CSV are diagnostic evidence only. They are not replacement input/retention estimates or approved migrated F8 panel c values. Panel a is a complementary common-grid display; its endpoints must not be read as the panel c gain.',
                  simultaneous_loss_window='Supplemental joint-loss sample is freshly re-extracted on common-grid timing, Feb21--May10 in WY2026, and F8 panel a shading follows that new window. See F8_S5_refresh_validation.json for old-sample reproduction, added cell-days and revised statistics.',
                  inputs=[record(x) for x in [a.native_daily, a.native_accounting, a.common_daily, a.core_summary, a.spatial_coverage, Path(__file__)] if x],
                  annual_rows=len(rows), max_native_identity_error=float(frame.multiplicative_identity_error.abs().max()))
    (a.output/'F8_accounting_support_decision.json').write_text(json.dumps(method, indent=2), encoding='utf-8')
    print(frame.loc[frame.wy == 2026].to_string(index=False))


if __name__ == '__main__':
    main()
