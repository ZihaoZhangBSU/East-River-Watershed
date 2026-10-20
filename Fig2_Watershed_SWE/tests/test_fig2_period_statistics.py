"""Synthetic checks distinguishing complete record periods from seasonal subsets."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Fig2_Watershed_SWE import model_statistics


class RecordPeriodStatisticsTests(unittest.TestCase):
    def test_historical_summary_includes_july_september_and_zero_days(self):
        dates = pd.date_range("2017-10-01", "2018-09-30")
        summer = dates >= "2018-07-01"
        table = pd.DataFrame({"date": dates, "wy": 2018,
                              "swe_i_mm": np.where(summer, 0., 2.),
                              "swe_n_mm": np.where(summer, 2., 1.)})
        result = model_statistics(table, pd.DataFrame({"wy": []})).iloc[0]
        self.assertEqual(result.n_dates, 365)
        self.assertEqual(result.end_date, "2018-09-30")
        self.assertAlmostEqual(result.mean_swe_isnobal_mm, 546 / 365)
        self.assertAlmostEqual(result.mean_swe_noahmp_mm, 457 / 365)
        self.assertAlmostEqual(result.relative_difference_noah_vs_isnobal_pct, 100 * (457 / 546 - 1))
        self.assertNotAlmostEqual(result.relative_difference_noah_vs_isnobal_pct, -50.)

    def test_wy2026_remains_the_available_october_june_record(self):
        dates = pd.date_range("2025-10-01", "2026-06-30")
        table = pd.DataFrame({"date": dates, "wy": 2026,
                              "swe_i_mm": np.arange(273, dtype=float),
                              "swe_n_mm": np.arange(273, dtype=float) / 2})
        result = model_statistics(table, pd.DataFrame({"wy": []})).iloc[0]
        self.assertEqual(result.n_dates, 273)
        self.assertEqual(result.end_date, "2026-06-30")
        self.assertAlmostEqual(result.relative_difference_noah_vs_isnobal_pct, -50.)
        self.assertAlmostEqual(result.mean_bias_noah_minus_isnobal_mm, -68.)

    def test_zero_reference_mean_has_no_relative_percentage(self):
        table = pd.DataFrame({"date": pd.date_range("2025-10-01", periods=2),
                              "wy": 2026, "swe_i_mm": [0., 0.], "swe_n_mm": [1., 2.]})
        with np.errstate(invalid="ignore", divide="ignore"):
            result = model_statistics(table, pd.DataFrame({"wy": []})).iloc[0]
        self.assertTrue(np.isnan(result.relative_difference_noah_vs_isnobal_pct))


if __name__ == "__main__":
    unittest.main()
