"""Meaningful synthetic checks for calendar support and diagnostic definitions."""
from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

import Station_Precip_SWE_Bias_and_Storage_Density as diagnostic


def complete_year(wy=2020):
    dates = diagnostic.year_dates(wy)
    return pd.DataFrame({"date": dates, "swe_mm": np.zeros(len(dates)),
                         "precip_mm": np.ones(len(dates)),
                         "source_precip_is_gapfilled": np.zeros(len(dates), dtype=bool)})


class DiagnosticDefinitions(unittest.TestCase):
    def test_calendar_keeps_leap_day_and_partial_year(self):
        self.assertEqual(len(diagnostic.year_dates(2020)), 366)
        self.assertIn(pd.Timestamp("2020-02-29"), diagnostic.year_dates(2020))
        self.assertEqual(diagnostic.year_dates(2026)[-1], pd.Timestamp("2026-06-30"))
        self.assertEqual(len(diagnostic.year_dates(2026)), 273)

    def test_bias_keeps_zero_days_and_uses_ratio_of_means(self):
        result = diagnostic.paired_bias(np.array([0, 1, 3, 0]), np.array([2, 2, 2, 2]))
        self.assertEqual(result["n_paired_days"], 4)
        self.assertEqual(result["mean_bias"], 1)
        self.assertEqual(result["relative_mean_bias_pct"], 100)

    def test_zero_reference_bias_is_explicitly_undefined(self):
        result = diagnostic.paired_bias(np.zeros(3), np.ones(3))
        self.assertTrue(np.isnan(result["relative_mean_bias_pct"]))
        self.assertEqual(result["relative_bias_status"], "zero_reference_mean")

    def test_bias_rejects_missing_or_unpaired_days(self):
        for reference, model in [(np.array([0, np.nan]), np.array([1, 2])),
                                 (np.array([0, 1]), np.array([1])),
                                 (np.array([]), np.array([]))]:
            with self.assertRaises(ValueError):
                diagnostic.paired_bias(reference, model)

    def test_density_exact_threshold_and_daily_mean(self):
        result = diagnostic.mean_density(np.array([10, 40, 999, 0, 50]),
                                          np.array([100, 200, 99.999, 0, np.nan]))
        self.assertEqual(result["density_n_days"], 2)
        self.assertEqual(result["density_mean_kg_m3"], 150)
        self.assertEqual(result["density_below_100mm_days"], 2)
        self.assertEqual(result["density_nonfinite_pair_days"], 1)

    def test_density_has_no_upper_cap_and_rejects_negative_values(self):
        self.assertEqual(diagnostic.mean_density(np.array([200]), np.array([100]))["density_mean_kg_m3"], 2000)
        with self.assertRaises(ValueError):
            diagnostic.mean_density(np.array([-1, 100]), np.array([50, 500]))

    def test_density_with_no_support_is_explicit(self):
        result = diagnostic.mean_density(np.array([0, 1]), np.array([0, 99]))
        self.assertTrue(np.isnan(result["density_mean_kg_m3"]))
        self.assertEqual(result["density_status"], "no_depth_support")

    def test_peak_first_exact_tie_and_excludes_july(self):
        d = complete_year()
        d.loc[d.date.isin(pd.to_datetime(["2020-02-29", "2020-03-01"])), "swe_mm"] = 100
        d.loc[d.date.eq(pd.Timestamp("2020-07-01")), "swe_mm"] = 200
        result = diagnostic.own_peak_ratio(d, 2020)
        self.assertEqual(result["peak_date"], "2020-02-29")
        self.assertEqual(result["peak_swe_mm"], 100)
        self.assertEqual(result["precip_to_peak_n_days"], 152)
        self.assertEqual(result["precip_to_own_peak_mm"], 152)
        self.assertAlmostEqual(result["peak_swe_precip_ratio"], 100/152)

    def test_sources_use_their_own_peak_precipitation_window(self):
        earlier, later = complete_year(), complete_year()
        earlier.loc[earlier.date.eq(pd.Timestamp("2020-02-01")), "swe_mm"] = 100
        later.loc[later.date.eq(pd.Timestamp("2020-04-01")), "swe_mm"] = 100
        a, b = diagnostic.own_peak_ratio(earlier, 2020), diagnostic.own_peak_ratio(later, 2020)
        self.assertLess(a["precip_to_peak_n_days"], b["precip_to_peak_n_days"])
        self.assertGreater(a["peak_swe_precip_ratio"], b["peak_swe_precip_ratio"])

    def test_peak_rejects_missing_precipitation_or_calendar_day(self):
        d = complete_year()
        d.loc[d.date.eq(pd.Timestamp("2020-02-29")), "precip_mm"] = np.nan
        with self.assertRaises(ValueError):
            diagnostic.own_peak_ratio(d, 2020)
        with self.assertRaises(ValueError):
            diagnostic.own_peak_ratio(complete_year().iloc[1:], 2020)

    def test_zero_precipitation_denominator_is_explicit(self):
        d = complete_year()
        d["precip_mm"] = 0
        result = diagnostic.own_peak_ratio(d, 2020)
        self.assertTrue(np.isnan(result["peak_swe_precip_ratio"]))
        self.assertEqual(result["peak_ratio_status"], "zero_precipitation")

    def test_gapfill_count_includes_only_dates_through_peak(self):
        d = complete_year()
        d.loc[d.date.eq(pd.Timestamp("2020-03-01")), "swe_mm"] = 10
        d.loc[d.date.isin(pd.to_datetime(["2020-02-29", "2020-04-01"])), "source_precip_is_gapfilled"] = True
        self.assertEqual(diagnostic.own_peak_ratio(d, 2020)["source_precip_gapfill_days_to_peak"], 1)


if __name__ == "__main__":
    unittest.main()
