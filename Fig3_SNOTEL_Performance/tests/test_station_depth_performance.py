"""Synthetic checks for station snow-depth scoring and figure semantics."""
import copy
import csv
import math
import os
import sys
import tempfile
import unittest
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import matplotlib.pyplot as plt
SCRIPT_DIR = Path(__file__).resolve().parent
if not (SCRIPT_DIR / 'Supp_Station_Depth_Performance.py').is_file():
    SCRIPT_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))
import Supp_Station_Depth_Performance as depth


def synthetic_sources():
    sources = {source: {} for source in depth.SOURCE_COLUMNS}
    for station in depth.STATIONS:
        for wy in depth.YEARS:
            for source in sources:
                values = {day: 0. for day in depth.amount_days(wy)}
                if source == 'SNOTEL':
                    values[date(wy, 3, 10)] = values[date(wy, 3, 11)] = 100.
                else:
                    values[date(wy, 3, 9 if source == 'iSnobal' else 13)] = 100. if source == 'iSnobal' else 120.
                    if wy != 2026:
                        values[date(wy, 7, 4)] = 500.  # Amount only; outside peak search.
                sources[source][station, wy] = values
    return sources


class StationDepthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources = synthetic_sources()
        cls.rows, cls.audit = depth.recompute(cls.sources)

    def test_zero_days_enter_the_amount_denominator(self):
        mb, rmse = depth.amount_metrics([0., 100., 0.], [0., 120., 0.])
        self.assertAlmostEqual(mb, 20 / 3)
        self.assertAlmostEqual(rmse, math.sqrt(400 / 3))
        with self.assertRaises(ValueError):
            depth.amount_metrics([0.], [float('nan')])
        with self.assertRaises(ValueError):
            depth.amount_metrics([0., 1.], [0.])

    def test_relative_font_is_normalized_before_registry_opens_it(self):
        supplied = Path('long-relative-folder/../fonts/Plot.ttf')
        actual = os.path.normpath(os.path.abspath(supplied))
        registry = SimpleNamespace(ttflist=[], addfont=mock.Mock(),
                                   _findfont_cached=SimpleNamespace(cache_clear=mock.Mock()))
        properties = mock.Mock()
        properties.get_name.return_value = 'Test Plot Font'
        with mock.patch.object(depth.font_manager, 'fontManager', registry), \
                mock.patch.object(depth.font_manager, 'FontProperties', return_value=properties), \
                mock.patch.object(depth.font_manager, 'findfont', return_value=actual):
            self.assertEqual(depth.configure_font(supplied), 'Test Plot Font')
            registry.addfont.assert_called_once_with(actual)

    def test_a_different_resolved_font_is_rejected(self):
        registry = SimpleNamespace(ttflist=[], addfont=mock.Mock(),
                                   _findfont_cached=SimpleNamespace(cache_clear=mock.Mock()))
        properties = mock.Mock()
        properties.get_name.return_value = 'Test Plot Font'
        with mock.patch.object(depth.font_manager, 'fontManager', registry), \
                mock.patch.object(depth.font_manager, 'FontProperties', return_value=properties), \
                mock.patch.object(depth.font_manager, 'findfont', return_value='different-font.ttf'):
            with self.assertRaisesRegex(ValueError, 'does not resolve'):
                depth.configure_font(Path('fonts/Plot.ttf'))

    def test_exact_tie_first_date_and_july_exclusion(self):
        values = {date(2019, 7, 1): 900., date(2019, 3, 11): 100.,
                  date(2019, 3, 9): 100. - 1e-12, date(2019, 3, 10): 100.}
        peak, maximum, ties = depth.first_peak(values, 2019)
        self.assertEqual((peak, maximum), (date(2019, 3, 10), 100.))
        self.assertEqual(ties, [date(2019, 3, 10), date(2019, 3, 11)])

    def test_depth_peaks_use_source_dates_and_model_minus_observation(self):
        for row in self.rows:
            self.assertEqual(row['observed_peak_date'], f"{row['wy']}-03-10")
            self.assertEqual(row['peak_error_days'], -1 if row['model'] == 'iSnobal' else 3)
            self.assertEqual(row['observed_peak_tie_count'], 2)
        self.assertEqual(len(self.audit['source_peak_events']), 81)

    def test_full_historical_amounts_but_partial_2026_and_bounded_peaks(self):
        rows = {(r['station_id'], r['wy'], r['model']): r for r in self.rows}
        for wy in depth.YEARS:
            row = rows[380, wy, 'Noah-MP']
            expected_n = 273 if wy == 2026 else 366 if wy in (2020, 2024) else 365
            self.assertEqual(row['amount_n_days'], expected_n)
            self.assertAlmostEqual(row['mb_mm'], (-80 if wy == 2026 else 420) / expected_n)
            self.assertEqual(row['peak_window_end'], f'{wy}-06-30')
            self.assertEqual(row['model_peak_date'], f'{wy}-03-13')
        self.assertEqual(len(self.rows), 54)
        self.assertFalse(any('sdd' in field.lower() or 'disappearance' in field.lower()
                             for field in depth.METRIC_FIELDS))

    def test_gap_in_any_source_is_rejected_without_filling_or_pair_filtering(self):
        changed = copy.deepcopy(self.sources)
        del changed['Noah-MP'][380, 2018][date(2018, 1, 1)]
        with self.assertRaisesRegex(ValueError, 'identical complete'):
            depth.recompute(changed)

    def test_saved_metric_wrong_period_count_or_peak_sign_is_rejected(self):
        for field, value in [('amount_window_end', '2018-06-30'), ('amount_n_days', 273),
                             ('peak_error_days', 1), ('observed_peak_date', '2018-07-01')]:
            changed = copy.deepcopy(self.rows)
            changed[0][field] = value
            with self.assertRaises(ValueError):
                depth.validate_metrics(changed)
        with self.assertRaisesRegex(ValueError, '54'):
            depth.validate_metrics(self.rows[:-1] + [self.rows[0]])

    def test_loader_uses_depth_and_rejects_duplicate_missing_and_nonfinite_values(self):
        fields = ['date', 'water_year', 'station_id', 'station_name', 'snow_depth_mm', 'swe_mm']
        rows = [{'date': day.isoformat(), 'water_year': wy, 'station_id': station,
                 'station_name': depth.STATIONS[station], 'snow_depth_mm': value, 'swe_mm': 99999.}
                for (station, wy), daily in self.sources['SNOTEL'].items()
                for day, value in daily.items()]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'daily.csv'
            def write(values):
                with path.open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(values)
            write(rows)
            loaded = depth.load_source(path, 'SNOTEL')
            self.assertEqual(loaded[380, 2018][date(2018, 3, 10)], 100.)
            write(rows + [rows[0]])
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                depth.load_source(path, 'SNOTEL')
            write(rows[1:])
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                depth.load_source(path, 'SNOTEL')
            for bad in ('nan', '-1', ''):
                changed = [dict(rows[0], snow_depth_mm=bad)] + rows[1:]
                write(changed)
                with self.assertRaises(ValueError):
                    depth.load_source(path, 'SNOTEL')

    def test_figure_has_nine_panels_162_markers_and_current_f3_model_styles(self):
        fig, audit = depth.make_figure(self.rows)
        try:
            self.assertEqual(len(fig.axes), 9)
            self.assertEqual(audit['marker_count'], 162)
            self.assertEqual(audit['connector_count'], 81)
            self.assertEqual({r['metric'] for r in audit['plotted_markers']}, set(depth.PLOT_FIELDS))
            self.assertEqual({r['marker'] for r in audit['plotted_markers'] if r['model'] == 'iSnobal'}, {'o'})
            self.assertEqual({r['marker'] for r in audit['plotted_markers'] if r['model'] == 'Noah-MP'}, {'s'})
            self.assertEqual(audit['minimum_font_pt'], 7)
            self.assertEqual(audit['font_hierarchy_pt'],
                             {'labels_titles_panel_letters': 8, 'ticks': 7, 'legend': 7})
            self.assertEqual({text.get_fontsize() for text in fig.legends[0].get_texts()}, {7})
            for axis in fig.axes:
                self.assertEqual(len(axis.patches), 1)  # WY2026 shading.
                self.assertEqual(sum(line.get_marker() in ('o', 's') for line in axis.lines), 18)
                self.assertEqual(sum(len(line.get_xdata()) == 2 and len(line.get_ydata()) == 2
                                     and line.get_ydata()[0] == line.get_ydata()[1]
                                     for line in axis.lines), 9)
        finally:
            plt.close(fig)


if __name__ == '__main__':
    unittest.main()
