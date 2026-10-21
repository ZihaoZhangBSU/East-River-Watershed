"""Synthetic input-contract tests; no research data or private paths required."""
import unittest

import numpy as np
import pandas as pd

from Fig4_SNOTEL_WY2026 import STATIONS, validate_inputs


def valid_inputs():
    dates = pd.date_range("2025-10-01", "2026-06-30")
    daily = pd.concat([
        pd.DataFrame({"date": dates, "station_id": station, "SNOTEL": 100.,
                      "iSnobal": 90., "Noah-MP": 80., "lo": 70., "hi": 110.})
        for station in STATIONS
    ], ignore_index=True)
    timing = pd.DataFrame([
        {"station_id": station, "wy": 2026, "source": source,
         "peak_date": "2026-03-01", "sdd_date": "2026-05-01"}
        for station in STATIONS for source in ["SNOTEL", "iSnobal", "Noah-MP"]
    ])
    return daily, timing


class Fig4InputTests(unittest.TestCase):
    def test_complete_three_station_support_is_accepted(self):
        daily, timing = valid_inputs()
        self.assertEqual(len(validate_inputs(daily, timing)), 9)

    def test_duplicate_date_is_rejected_at_unchanged_row_count(self):
        daily, timing = valid_inputs()
        daily.loc[1, "date"] = daily.loc[0, "date"]
        with self.assertRaises(ValueError):
            validate_inputs(daily, timing)

    def test_missing_calendar_day_is_rejected(self):
        daily, timing = valid_inputs()
        with self.assertRaises(ValueError):
            validate_inputs(daily.drop(index=1), timing)

    def test_nonfinite_curve_or_range_value_is_rejected(self):
        for column in ["SNOTEL", "iSnobal", "Noah-MP", "lo", "hi"]:
            with self.subTest(column=column):
                daily, timing = valid_inputs()
                daily.loc[1, column] = np.nan
                with self.assertRaises(ValueError):
                    validate_inputs(daily, timing)

    def test_inverted_range_is_rejected(self):
        daily, timing = valid_inputs()
        daily.loc[1, ["lo", "hi"]] = [120., 100.]
        with self.assertRaises(ValueError):
            validate_inputs(daily, timing)

    def test_curve_outside_range_is_rejected(self):
        daily, timing = valid_inputs()
        daily.loc[1, "iSnobal"] = 120.
        with self.assertRaises(ValueError):
            validate_inputs(daily, timing)

    def test_incomplete_timing_schema_is_rejected(self):
        daily, timing = valid_inputs()
        with self.assertRaises(ValueError):
            validate_inputs(daily, timing.drop(columns="peak_date"))

    def test_missing_timing_key_is_rejected(self):
        daily, timing = valid_inputs()
        with self.assertRaises(ValueError):
            validate_inputs(daily, timing.drop(index=0))

    def test_duplicate_timing_key_is_rejected_at_unchanged_row_count(self):
        daily, timing = valid_inputs()
        timing.loc[1] = timing.loc[0]
        with self.assertRaises(ValueError):
            validate_inputs(daily, timing)

    def test_disappearance_on_or_before_own_peak_is_rejected(self):
        for invalid_date in ["2026-03-01", "2026-02-28"]:
            with self.subTest(sdd_date=invalid_date):
                daily, timing = valid_inputs()
                timing.loc[0, "sdd_date"] = invalid_date
                with self.assertRaises(ValueError):
                    validate_inputs(daily, timing)


if __name__ == "__main__":
    unittest.main()
