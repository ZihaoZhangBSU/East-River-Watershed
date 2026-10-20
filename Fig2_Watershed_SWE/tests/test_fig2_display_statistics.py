"""Independent synthetic definitions for Figure 2's daily intermodel annotations."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Fig2_Watershed_SWE import paired_model_statistics, statistics_annotation


class DisplayStatisticsTests(unittest.TestCase):
    def test_multiplier_with_finite_zero_day(self):
        result = paired_model_statistics([0, 1, 2, 3], [0, 2, 4, 6])
        self.assertAlmostEqual(result["pearson_r"], 1.)
        self.assertEqual(result["mean_bias_noah_minus_isnobal_mm"], 1.5)
        self.assertEqual(result["relative_difference_noah_vs_isnobal_pct"], 100.)
        self.assertEqual(result["mean_swe_isnobal_mm"], 1.5)

    def test_sign_and_reference_are_not_reversed(self):
        result = paired_model_statistics([1, 2, 3], [2, 3, 4])
        self.assertEqual(result["mean_bias_noah_minus_isnobal_mm"], 1.)
        self.assertEqual(result["relative_difference_noah_vs_isnobal_pct"], 50.)
        reverse = paired_model_statistics([2, 3, 4], [1, 2, 3])
        self.assertEqual(reverse["mean_bias_noah_minus_isnobal_mm"], -1.)
        self.assertAlmostEqual(reverse["relative_difference_noah_vs_isnobal_pct"], -100 / 3)

    def test_negative_correlation_and_equal_means(self):
        result = paired_model_statistics([0, 1, 2, 3], [3, 2, 1, 0])
        self.assertAlmostEqual(result["pearson_r"], -1.)
        self.assertEqual(result["mean_bias_noah_minus_isnobal_mm"], 0.)
        self.assertEqual(result["relative_difference_noah_vs_isnobal_pct"], 0.)

    def test_constant_nonzero_reference_has_undefined_r_but_defined_relative_mb(self):
        result = paired_model_statistics([1, 1, 1], [1, 2, 3])
        self.assertTrue(np.isnan(result["pearson_r"]))
        self.assertEqual(result["relative_difference_noah_vs_isnobal_pct"], 100.)

    def test_either_constant_series_makes_r_undefined(self):
        for i, n in [([1, 2, 3], [1, 1, 1]), ([1], [3])]:
            self.assertTrue(np.isnan(paired_model_statistics(i, n)["pearson_r"]))

    def test_all_zero_does_not_invent_agreement_or_relative_bias(self):
        result = paired_model_statistics([0, 0], [0, 0])
        self.assertTrue(np.isnan(result["pearson_r"]))
        self.assertEqual(result["mean_bias_noah_minus_isnobal_mm"], 0.)
        self.assertTrue(np.isnan(result["relative_difference_noah_vs_isnobal_pct"]))
        annotation = statistics_annotation(pd.Series(result))
        self.assertEqual(annotation, "r = NA\nMB = +0.0 mm\nrMB = NA")

    def test_zero_reference_with_positive_other_model_stays_undefined(self):
        result = paired_model_statistics([0, 0], [1, 3])
        self.assertEqual(result["mean_bias_noah_minus_isnobal_mm"], 2.)
        self.assertTrue(np.isnan(result["relative_difference_noah_vs_isnobal_pct"]))

    def test_missing_or_infinite_pairs_cannot_change_statistics_support(self):
        for bad in [np.nan, np.inf, -np.inf]:
            with self.assertRaisesRegex(ValueError, "finite pairs"):
                paired_model_statistics([1, bad], [2, 3])
            with self.assertRaisesRegex(ValueError, "finite pairs"):
                paired_model_statistics([1, 2], [2, bad])

    def test_empty_unequal_or_multidimensional_series_are_rejected(self):
        for i, n in [([], []), ([1], [1, 2]), ([[1, 2]], [[1, 2]])]:
            with self.assertRaisesRegex(ValueError, "1D paired"):
                paired_model_statistics(i, n)

    def test_annotation_uses_declared_rounding_and_units(self):
        result = paired_model_statistics([0, 1, 2, 3], [0, 2, 4, 6])
        self.assertEqual(statistics_annotation(pd.Series(result)),
                         "r = 1.000\nMB = +1.5 mm\nrMB = +100.0%")


if __name__ == "__main__":
    unittest.main()
