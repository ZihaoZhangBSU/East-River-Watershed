"""Synthetic safety checks for the presentation-only elevation helper."""
import copy
import csv
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import Supp_ASO_SWE_Elevation as figure


def synthetic_rows():
    rows = []
    for date_index, date in enumerate(figure.FLIGHT_DATES):
        for model_index, model in enumerate(figure.MODELS):
            for band_index, band in enumerate(figure.ELEVATION_BANDS):
                rows.append({'date': date, 'variable': 'SWE', 'model': model,
                             'elevation_band': band, 'n_cells': 20 + band_index,
                             'mb_mm': -date_index - 10 * model_index - band_index,
                             'rmse_mm': 10 + date_index + model_index + band_index,
                             'r': (date_index + band_index) / 20})
    return rows


class ElevationFigureTests(unittest.TestCase):
    def font_registry_mocks(self, actual):
        registry = SimpleNamespace(ttflist=[], addfont=mock.Mock(),
                                   _findfont_cached=SimpleNamespace(cache_clear=mock.Mock()))
        properties = mock.Mock()
        properties.get_name.return_value = 'Test Plot Font'
        return (mock.patch.object(figure.font_manager, 'fontManager', registry),
                mock.patch.object(figure.font_manager, 'FontProperties', return_value=properties),
                mock.patch.object(figure.font_manager, 'findfont', return_value=str(actual)))

    def test_relative_font_dotdot_matches_normalized_registry_path(self):
        supplied = Path('relative-fonts/../plot-fonts/Plot.ttf')
        actual = os.path.normpath(os.path.abspath(supplied))
        registry, properties, findfont = self.font_registry_mocks(actual)
        # Reproduce a sandbox resolve() that retains the supplied spelling.
        with registry, properties, findfont, mock.patch.object(Path, 'resolve', lambda self: self):
            self.assertEqual(figure.configure_font(supplied), 'Test Plot Font')

    def test_different_font_file_is_still_rejected(self):
        supplied = Path('relative-fonts/../plot-fonts/Plot.ttf')
        actual = os.path.normpath(os.path.abspath('plot-fonts/Different.ttf'))
        registry, properties, findfont = self.font_registry_mocks(actual)
        with registry, properties, findfont:
            with self.assertRaisesRegex(ValueError, 'did not resolve'):
                figure.configure_font(supplied)

    def test_date_model_band_mapping_survives_reversed_input(self):
        arrays = figure.score_arrays(list(reversed(synthetic_rows())))
        for di in range(14):
            for mi, model in enumerate(figure.MODELS):
                for bi in range(3):
                    self.assertEqual(arrays[model, 'mb_mm'][di, bi], -di - 10 * mi - bi)
        self.assertEqual(arrays['iSnobal', 'r'].shape, (14, 3))

    def test_only_swe_named_bands_are_selected(self):
        rows = synthetic_rows()
        rows.extend([dict(rows[0], variable='Depth'), dict(rows[0], elevation_band='All')])
        self.assertEqual(len(figure.select_band_rows(rows)), 84)

    def test_duplicate_or_missing_keys_are_rejected(self):
        rows = synthetic_rows()
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            figure.select_band_rows(rows + [copy.deepcopy(rows[0])])
        with self.assertRaisesRegex(ValueError, '84'):
            figure.select_band_rows(rows[:-1])

    def test_undefined_r_remains_missing_and_finite_zero_stays_zero(self):
        rows = synthetic_rows()
        rows[3]['r'] = None
        arrays = figure.score_arrays(rows)
        self.assertTrue(np.isnan(arrays['Noah-MP', 'r'][0, 0]))
        self.assertEqual(arrays['iSnobal', 'r'][0, 0], 0.)
        self.assertFalse(np.ma.masked_invalid(arrays['iSnobal', 'r']).mask[0, 0])
        self.assertTrue(np.ma.masked_invalid(arrays['Noah-MP', 'r']).mask[0, 0])

    def test_shared_scales_are_signed_and_prevent_silent_clipping(self):
        arrays = figure.score_arrays(synthetic_rows())
        self.assertEqual(figure.metric_scales(arrays),
                         {'mb_mm': (-300., 300.), 'rmse_mm': (0., 400.), 'r': (-1., 1.)})
        for metric, outside in [('mb_mm', -301), ('rmse_mm', 401), ('r', 1.01)]:
            changed = {key: value.copy() for key, value in arrays.items()}
            changed['Noah-MP', metric][13, 2] = outside
            with self.assertRaisesRegex(ValueError, 'exceeds'):
                figure.metric_scales(changed)

    def test_export_preserves_scalar_strings_and_support(self):
        fields = ['date', 'variable', 'model', 'elevation_band', 'n_cells', 'r']
        rows = [{'date': '2018-05-24', 'variable': 'SWE', 'model': 'Noah-MP',
                 'elevation_band': 'Lower', 'n_cells': '187', 'r': ''},
                {'date': '2018-05-24', 'variable': 'SWE', 'model': 'iSnobal',
                 'elevation_band': 'Upper', 'n_cells': '225', 'r': '0.9402659013645762'},
                {'date': '2018-05-24', 'variable': 'Depth', 'model': 'iSnobal',
                 'elevation_band': 'Upper', 'n_cells': '222', 'r': '0.9'}]
        with tempfile.TemporaryDirectory() as directory:
            source, output = [Path(directory) / name for name in ('in.csv', 'out.csv')]
            with source.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            self.assertEqual(figure.export_exact_band_rows(source, output), 2)
            with output.open(newline='') as stream:
                self.assertEqual(list(csv.DictReader(stream)), rows[:2])

    def test_json_arrays_use_null_instead_of_nan(self):
        self.assertEqual(figure.json_array(np.array([[0., np.nan, -0.2]])),
                         [[0., None, -0.2]])


if __name__ == '__main__':
    unittest.main()
