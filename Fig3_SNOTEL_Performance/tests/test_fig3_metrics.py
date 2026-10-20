"""Scientific edge cases and strict source validation; no research data needed."""
import csv
import math
import tempfile
import unittest
from datetime import date, timedelta
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import recompute_fig3_metrics as metrics


class MetricDefinitionTests(unittest.TestCase):
    def test_amount_sign_and_zero_inclusion(self):
        bias, rmse = metrics.amount_metrics([0.0, 10.0], [2.0, 6.0])
        self.assertEqual(bias, -1.0)
        self.assertAlmostEqual(rmse, math.sqrt(10))

    def test_amount_rejects_unpaired_or_nonfinite(self):
        for obs, model in [([], []), ([0], [0, 1]), ([0], [float("nan")])]:
            with self.assertRaises(ValueError):
                metrics.amount_metrics(obs, model)

    def test_peak_first_exact_tie_and_october_june_bounds(self):
        daily = {date(2019, 9, 30): 900, date(2019, 10, 1): 2,
                 date(2020, 2, 29): 10, date(2020, 3, 1): 10,
                 date(2020, 7, 1): 999}
        peak, value, ties = metrics.first_peak(daily, 2020)
        self.assertEqual(peak, date(2020, 2, 29))
        self.assertEqual(value, 10)
        self.assertEqual(ties, [date(2020, 2, 29), date(2020, 3, 1)])

    def test_all_zero_series_retains_own_peak_and_starts_after_it(self):
        first = date(2025, 10, 1)
        daily = {first + timedelta(days=i): 0.0 for i in range(8)}
        peak, value, _ = metrics.first_peak(daily, 2026)
        self.assertEqual((peak, value), (first, 0))
        self.assertEqual(metrics.snow_disappearance(daily, peak), first + timedelta(days=1))

    def test_exact_10_breaks_run_and_sdd_is_first_not_seventh_day(self):
        peak = date(2026, 5, 1)
        daily = {peak + timedelta(days=i): 0 for i in range(11)}
        daily[peak] = 100
        daily[peak + timedelta(days=3)] = 10
        self.assertEqual(metrics.snow_disappearance(daily, peak), peak + timedelta(days=4))

    def test_missing_date_and_nonfinite_break_seven_calendar_days(self):
        peak = date(2026, 5, 1)
        for bad in (None, float("nan"), float("inf")):
            daily = {peak + timedelta(days=i): 0 for i in range(12)}
            daily[peak] = 100
            if bad is None:
                del daily[peak + timedelta(days=4)]
            else:
                daily[peak + timedelta(days=4)] = bad
            self.assertEqual(metrics.snow_disappearance(daily, peak), peak + timedelta(days=5))

    def test_six_day_endpoint_is_censored(self):
        peak = date(2026, 6, 24)
        daily = {peak + timedelta(days=i): 0 for i in range(7)}
        daily[peak] = 100
        self.assertIsNone(metrics.snow_disappearance(daily, peak))

    def test_historical_sdd_can_be_confirmed_after_june(self):
        peak = date(2024, 6, 29)
        daily = {peak + timedelta(days=i): 0 for i in range(9)}
        daily[peak] = 100
        daily[peak + timedelta(days=1)] = 20
        self.assertEqual(metrics.snow_disappearance(daily, peak), date(2024, 7, 1))

    def test_post_own_peak_not_observed_peak(self):
        first = date(2026, 5, 1)
        daily = {first + timedelta(days=i): 0 for i in range(20)}
        own_peak = first + timedelta(days=10)
        daily[own_peak] = 100
        self.assertEqual(metrics.snow_disappearance(daily, own_peak), own_peak + timedelta(days=1))

    def test_leap_scoring_calendars(self):
        for wy, expected in [(2019, 273), (2020, 274), (2024, 274), (2026, 273)]:
            self.assertEqual(len(list(metrics.calendar_days(date(wy - 1, 10, 1), date(wy, 6, 30)))), expected)


class SourceValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "daily.csv"
        self.rows = [{"date": day.isoformat(), "water_year": wy,
                      "station_id": station, "station_name": name, "swe_mm": 0.0}
                     for station, name in metrics.STATIONS.items() for wy in metrics.YEARS
                     for day in metrics.calendar_days(date(wy - 1, 10, 1), metrics.available_end(wy))]

    def tearDown(self):
        self.temp.cleanup()

    def write(self, rows):
        with self.path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["date", "water_year", "station_id", "station_name", "swe_mm"])
            writer.writeheader()
            writer.writerows(rows)

    def test_complete_source_passes(self):
        self.write(self.rows)
        loaded = metrics.load_source(self.path, "SNOTEL")
        self.assertEqual(len(loaded), 27)
        self.assertEqual(sum(map(len, loaded.values())), 9585)

    def test_duplicate_replacing_missing_day_fails_even_with_expected_row_count(self):
        rows = [row.copy() for row in self.rows]
        rows[1] = rows[0].copy()
        self.write(rows)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            metrics.load_source(self.path, "SNOTEL")

    def test_missing_calendar_day_fails(self):
        self.write(self.rows[1:])
        with self.assertRaisesRegex(ValueError, "incomplete calendar"):
            metrics.load_source(self.path, "SNOTEL")

    def test_invalid_numeric_and_metadata_fail(self):
        variants = [("swe_mm", "nan"), ("swe_mm", "inf"), ("swe_mm", -1),
                    ("station_id", 999), ("water_year", 2019),
                    ("station_name", "Wrong station"), ("date", "2017-10-01T23:00:00")]
        for field, value in variants:
            with self.subTest(field=field, value=value):
                rows = [row.copy() for row in self.rows]
                rows[0][field] = value
                self.write(rows)
                with self.assertRaises(ValueError):
                    metrics.load_source(self.path, "SNOTEL")

    def test_recompute_censoring_and_zero_support(self):
        sources = {}
        for source in metrics.SOURCE_COLUMNS:
            sources[source] = {}
            for station in metrics.STATIONS:
                for wy in metrics.YEARS:
                    days = list(metrics.calendar_days(date(wy - 1, 10, 1), metrics.available_end(wy)))
                    sources[source][(station, wy)] = {day: 0.0 for day in days}
        # Observed unresolved, both models resolved: no invented SDD error.
        sources["SNOTEL"][(380, 2026)] = {day: 50.0 for day in sources["SNOTEL"][(380, 2026)]}
        rows, audit = metrics.recompute(sources)
        self.assertEqual(len(rows), 54)
        self.assertEqual(len(audit["source_events"]), 81)
        selected = [row for row in rows if row["station_id"] == 380 and row["wy"] == 2026]
        self.assertEqual([row["sdd_error_days"] for row in selected], ["", ""])
        self.assertTrue(all(not row["sdd_censored"] for row in selected))
        observed_event = next(row for row in audit["source_events"] if row["station_id"] == 380
                              and row["wy"] == 2026 and row["source"] == "SNOTEL")
        self.assertTrue(observed_event["sdd_censored"])
        self.assertEqual(sum(row["paired_scoring_days"] for row in audit["paired_support"]), 7377)
        self.assertEqual(sum(row["observed_zero_swe_days"] for row in audit["paired_support"]), 7104)


if __name__ == "__main__":
    unittest.main()
