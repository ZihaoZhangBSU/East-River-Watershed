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


def synthetic_rows(variable='SWE', offset=0):
    rows = []
    for date_index, date in enumerate(figure.FLIGHT_DATES):
        for model_index, model in enumerate(figure.MODELS):
            for band_index, band in enumerate(figure.ELEVATION_BANDS):
                rows.append({'date': date, 'variable': variable, 'model': model,
                             'elevation_band': band, 'n_cells': 20 + band_index + offset,
                             'mb_mm': -date_index - 10 * model_index - band_index - offset,
                             'rmse_mm': 10 + date_index + model_index + band_index + offset,
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
        arrays = figure.score_arrays(list(reversed(synthetic_rows() + synthetic_rows('Depth', 100))))
        for di in range(14):
            for mi, model in enumerate(figure.MODELS):
                for bi in range(3):
                    self.assertEqual(arrays[model, 'mb_mm'][di, bi], -di - 10 * mi - bi)
        self.assertEqual(arrays['iSnobal', 'r'].shape, (14, 3))

    def test_only_swe_named_bands_are_selected(self):
        rows = synthetic_rows()
        rows.extend([dict(rows[0], variable='Depth'), dict(rows[0], elevation_band='All')])
        self.assertEqual(len(figure.select_band_rows(rows)), 84)

    def test_depth_uses_its_own_values_and_support(self):
        rows = synthetic_rows() + synthetic_rows('Depth', offset=100)
        rows[85]['r'] = None  # A depth-only NA must not be replaced by SWE support.
        selected = figure.select_band_rows(rows, 'Depth')
        self.assertEqual(len(selected), 84)
        self.assertEqual({r['variable'] for r in selected}, {'Depth'})
        self.assertEqual({r['n_cells'] for r in selected}, {120, 121, 122})
        arrays = figure.score_arrays(rows, 'Depth')
        self.assertEqual(arrays['iSnobal', 'mb_mm'][0, 0], -100.)
        self.assertEqual(arrays['Noah-MP', 'rmse_mm'][13, 2], 126.)
        self.assertTrue(np.isnan(arrays['iSnobal', 'r'][0, 1]))
        self.assertEqual(figure.score_arrays(rows, 'SWE')['iSnobal', 'r'][0, 1], .05)
        with self.assertRaisesRegex(ValueError, '84 Depth'):
            figure.select_band_rows(rows[:-1], 'Depth')
        with self.assertRaisesRegex(ValueError, 'Unknown'):
            figure.select_band_rows(rows, 'both')

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
        for variable, mb_limit, rmse_limit in [('SWE', 300., 400.), ('Depth', 600., 800.)]:
            arrays = figure.score_arrays(synthetic_rows(variable), variable)
            scales = figure.metric_scales(arrays, variable)
            self.assertEqual(scales, {'mb_mm': (-mb_limit, mb_limit),
                                     'rmse_mm': (0., rmse_limit), 'r': (0., 1.)})
            r_norm = figure.colors.Normalize(*scales['r'])
            self.assertEqual(r_norm(0.), 0.)
            self.assertEqual(r_norm(1.), 1.)
            for metric, outside in [('mb_mm', -mb_limit - 1), ('mb_mm', mb_limit + 1),
                                    ('rmse_mm', rmse_limit + 1), ('r', 1.01), ('r', -.01)]:
                changed = {key: value.copy() for key, value in arrays.items()}
                changed['Noah-MP', metric][13, 2] = outside
                with self.assertRaisesRegex(ValueError, 'exceeds'):
                    figure.metric_scales(changed, variable)

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
            self.assertEqual(figure.export_exact_band_rows(source, output, 'Depth'), 1)
            with output.open(newline='') as stream:
                self.assertEqual(list(csv.DictReader(stream)), rows[2:])

    def test_json_arrays_use_null_instead_of_nan(self):
        self.assertEqual(figure.json_array(np.array([[0., np.nan, -0.2]])),
                         [[0., None, -0.2]])


if __name__ == '__main__':
    unittest.main()
