"""Synthetic Fig. 3 input and artist checks; contains no research data."""
from pathlib import Path
import itertools
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Fig3_SNOTEL_Performance import load_metrics, make_figure


def valid_table():
    names = {380: 'Butte', 737: 'Schofield Pass', 1141: 'Upper Taylor'}
    rows = []
    for sid, wy, model in itertools.product(names, range(2018, 2027), ['iSnobal', 'Noah-MP']):
        shift = 1 if model == 'iSnobal' else -1
        rows.append(dict(station_id=sid, station_name=names[sid], wy=wy, model=model,
                         mb=shift * (wy - 2017), rmse=wy - 2017 + 2,
                         peak_date=(pd.Timestamp(wy, 3, 1) + pd.Timedelta(days=shift)).strftime('%Y-%m-%d'),
                         sdd_date=(pd.Timestamp(wy, 5, 1) + pd.Timedelta(days=shift)).strftime('%Y-%m-%d'),
                         peak_error_days=shift, sdd_error_days=shift, sdd_censored=False,
                         amount_window_start=f'{wy-1}-10-01',
                         amount_window_end=f'{wy}-06-30' if wy == 2026 else f'{wy}-09-30',
                         amount_n_days=(pd.Timestamp(wy, 6, 30) - pd.Timestamp(wy-1, 10, 1)).days + 1 if wy == 2026
                                       else (pd.Timestamp(wy, 9, 30) - pd.Timestamp(wy-1, 10, 1)).days + 1,
                         available_record_end=f'{wy}-06-30' if wy == 2026 else f'{wy}-09-30'))
    return pd.DataFrame(rows).astype({'mb': float, 'rmse': float})


class Fig3InputTests(unittest.TestCase):
    def load(self, table):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'metrics.csv'
            table.to_csv(path, index=False)
            return load_metrics(path)

    def rejected(self, table):
        with self.assertRaises(ValueError):
            self.load(table)

    def test_complete_product_is_accepted(self):
        self.assertEqual(len(self.load(valid_table())), 54)

    def test_row_order_is_irrelevant(self):
        self.assertEqual(len(self.load(valid_table().sample(frac=1, random_state=4))), 54)

    def test_missing_model_row_is_rejected(self):
        self.rejected(valid_table().iloc[:-1])

    def test_duplicate_replacing_missing_row_is_rejected(self):
        table = valid_table()
        table.iloc[-1] = table.iloc[0]
        self.rejected(table)

    def test_extra_model_is_rejected(self):
        table = valid_table()
        table.loc[0, 'model'] = 'other'
        self.rejected(table)

    def test_unknown_station_is_rejected(self):
        table = valid_table()
        table.loc[0, 'station_id'] = 999
        self.rejected(table)

    def test_noninteger_water_year_is_rejected(self):
        table = valid_table().astype({'wy': float})
        table.loc[0, 'wy'] = 2018.5
        self.rejected(table)

    def test_missing_required_column_is_rejected(self):
        self.rejected(valid_table().drop(columns=['mb']))

    def test_nonfinite_amounts_are_rejected(self):
        for field, value in itertools.product(['mb', 'rmse'], [np.nan, np.inf, -np.inf]):
            with self.subTest(field=field, value=value):
                table = valid_table()
                table.loc[0, field] = value
                self.rejected(table)

    def test_historical_october_june_amount_metadata_is_rejected(self):
        table = valid_table()
        table.loc[0, 'amount_window_end'] = '2018-06-30'
        table.loc[0, 'amount_n_days'] = 273
        self.rejected(table)

    def test_2026_padded_full_year_metadata_is_rejected(self):
        table = valid_table()
        index = table.index[table.wy.eq(2026)][0]
        table.loc[index, 'amount_window_end'] = '2026-09-30'
        table.loc[index, 'amount_n_days'] = 365
        self.rejected(table)

    def test_amount_count_must_include_leap_day(self):
        table = valid_table()
        index = table.index[table.wy.eq(2020)][0]
        table.loc[index, 'amount_n_days'] = 365
        self.rejected(table)

    def test_negative_rmse_is_rejected(self):
        table = valid_table()
        table.loc[0, 'rmse'] = -1
        self.rejected(table)

    def test_rmse_below_absolute_bias_is_rejected(self):
        table = valid_table()
        table.loc[0, ['mb', 'rmse']] = [3, 2]
        self.rejected(table)

    def test_noninteger_timing_error_is_rejected(self):
        table = valid_table().astype({'peak_error_days': float})
        table.loc[0, 'peak_error_days'] = 0.5
        self.rejected(table)

    def test_nonfinite_peak_error_is_rejected(self):
        for value in [np.nan, np.inf, -np.inf]:
            with self.subTest(value=value):
                table = valid_table().astype({'peak_error_days': float})
                table.loc[0, 'peak_error_days'] = value
                self.rejected(table)

    def test_unknown_boolean_flag_is_rejected(self):
        table = valid_table().astype({'sdd_censored': object})
        table.loc[0, 'sdd_censored'] = 'sometimes'
        self.rejected(table)

    def test_invalid_calendar_date_is_rejected(self):
        table = valid_table()
        table.loc[0, 'peak_date'] = '2018-02-30'
        self.rejected(table)

    def test_peak_outside_october_june_is_rejected(self):
        table = valid_table()
        table.loc[0, 'peak_date'] = '2018-07-01'
        self.rejected(table)

    def test_incorrect_record_endpoint_is_rejected(self):
        table = valid_table()
        table.loc[0, 'available_record_end'] = '2018-06-30'
        self.rejected(table)

    def test_sdd_must_follow_own_peak(self):
        table = valid_table()
        table.loc[0, 'sdd_date'] = table.loc[0, 'peak_date']
        self.rejected(table)

    def test_unconfirmed_endpoint_sdd_is_rejected(self):
        table = valid_table()
        table.loc[0, 'sdd_date'] = '2018-09-25'
        self.rejected(table)

    def test_inconsistent_observed_peak_is_rejected(self):
        table = valid_table()
        table.loc[0, 'peak_error_days'] += 1
        self.rejected(table)

    def test_inconsistent_observed_sdd_is_rejected(self):
        table = valid_table()
        table.loc[0, 'sdd_error_days'] += 1
        self.rejected(table)

    def test_blank_sdd_with_no_censor_flag_is_rejected(self):
        table = valid_table()
        table.loc[0, ['sdd_date', 'sdd_error_days']] = [None, np.nan]
        self.rejected(table)

    def test_censored_model_event_keeps_error_missing(self):
        table = valid_table()
        table.loc[0, ['sdd_date', 'sdd_error_days', 'sdd_censored']] = [None, np.nan, True]
        loaded = self.load(table)
        row = loaded[(loaded.station_id == 380) & (loaded.wy == 2018) & (loaded.model == 'iSnobal')].iloc[0]
        self.assertTrue(pd.isna(row.sdd_error_days))

    def test_censored_error_is_not_zero(self):
        table = valid_table()
        table.loc[0, ['sdd_date', 'sdd_error_days', 'sdd_censored']] = [None, 0, True]
        self.rejected(table)

    def test_observed_censor_flag_allows_unavailable_error(self):
        table = valid_table()
        table['observed_sdd_censored'] = False
        mask = (table.station_id == 380) & (table.wy == 2018)
        table.loc[mask, 'observed_sdd_censored'] = True
        table.loc[mask, 'sdd_error_days'] = np.nan
        loaded = self.load(table)
        self.assertTrue(loaded.loc[(loaded.station_id == 380) & (loaded.wy == 2018), 'sdd_error_days'].isna().all())

    def test_observed_censor_flag_must_agree_between_models(self):
        table = valid_table()
        table['observed_sdd_censored'] = False
        table.loc[0, 'observed_sdd_censored'] = True
        table.loc[0, 'sdd_error_days'] = np.nan
        self.rejected(table)


class Fig3ArtistTests(unittest.TestCase):
    def test_all_zero_scores_have_usable_axes(self):
        table = valid_table()
        table.loc[:, ['mb', 'rmse', 'peak_error_days', 'sdd_error_days']] = 0
        # Identical model dates imply the same observed peak and SDD.
        table['peak_date'] = table.wy.map(lambda wy: f'{wy}-03-01')
        table['sdd_date'] = table.wy.map(lambda wy: f'{wy}-05-01')
        fig = make_figure(table)
        try:
            for ax in fig.axes:
                self.assertGreater(ax.get_xlim()[1], ax.get_xlim()[0])
                self.assertLessEqual(ax.get_xlim()[0], 0)
                self.assertGreaterEqual(ax.get_xlim()[1], 0)
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_all_censored_sdds_have_c_marks_and_finite_axes(self):
        table = valid_table()
        table['sdd_date'] = None
        table['sdd_error_days'] = np.nan
        table['sdd_censored'] = True
        table['observed_sdd_censored'] = True
        fig = make_figure(table)
        try:
            for panel in [3, 7, 11]:
                ax = fig.axes[panel]
                self.assertTrue(np.isfinite(ax.get_xlim()).all())
                self.assertGreater(ax.get_xlim()[1], ax.get_xlim()[0])
                self.assertEqual([text.get_text() for text in ax.texts].count('C'), 18)
                points = [line for line in ax.lines if line.get_marker() in ['o', 's']]
                self.assertEqual(points, [])
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_216_values_station_year_model_and_labels(self):
        table = valid_table()
        fig = make_figure(table)
        fields = ['mb', 'rmse', 'peak_error_days', 'sdd_error_days']
        stations = [(380, 'Butte'), (737, 'Schofield Pass'), (1141, 'Upper Taylor*')]
        count = 0
        try:
            self.assertEqual(len(fig.axes), 12)
            for row, (station_id, label) in enumerate(stations):
                self.assertIn(label, fig.axes[row * 4].get_ylabel())
                for col, field in enumerate(fields):
                    ax = fig.axes[row * 4 + col]
                    points = [line for line in ax.lines if line.get_marker() in ['o', 's']]
                    self.assertEqual(len(points), 18)
                    for year_index, wy in enumerate(range(2018, 2027)):
                        for model_index, (model, marker) in enumerate([('iSnobal', 'o'), ('Noah-MP', 's')]):
                            point = points[year_index * 2 + model_index]
                            expected = table[(table.station_id == station_id) & (table.wy == wy) & (table.model == model)].iloc[0][field]
                            self.assertEqual(float(point.get_xdata()[0]), expected)
                            self.assertEqual(float(point.get_ydata()[0]), year_index)
                            self.assertEqual(point.get_marker(), marker)
                            self.assertLessEqual(ax.get_xlim()[0], expected)
                            self.assertGreaterEqual(ax.get_xlim()[1], expected)
                            count += 1
            self.assertEqual(count, 216)
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_censored_sdd_is_marked_c_without_zero_point(self):
        table = valid_table()
        table.loc[0, ['sdd_date', 'sdd_error_days', 'sdd_censored']] = [None, np.nan, True]
        fig = make_figure(table)
        try:
            ax = fig.axes[3]
            points = [line for line in ax.lines if line.get_marker() in ['o', 's']]
            self.assertEqual(len(points), 17)
            self.assertEqual([text.get_text() for text in ax.texts].count('C'), 1)
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)


if __name__ == '__main__':
    unittest.main()
