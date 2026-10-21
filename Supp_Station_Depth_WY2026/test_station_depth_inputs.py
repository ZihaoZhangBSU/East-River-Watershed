"""Synthetic calendar and artist checks; no research data or private paths."""
import unittest

import numpy as np
import pandas as pd

import Supp_Station_Depth_WY2026 as depth


def valid_input():
    dates = pd.date_range("2025-10-01", "2026-06-30")
    return pd.concat([
        pd.DataFrame({"date": dates, "station_id": station, "SNOTEL": 100.,
                      "iSnobal": 90., "Noah-MP": 80., "lo": 70., "hi": 110.})
        for station in depth.STATIONS
    ], ignore_index=True)


class StationDepthTests(unittest.TestCase):
    def tearDown(self):
        depth.plt.close("all")

    def test_complete_support_is_accepted(self):
        self.assertIsNone(depth.validate_inputs(valid_input()))

    def test_missing_and_duplicate_dates_are_rejected(self):
        daily = valid_input()
        with self.assertRaises(ValueError):
            depth.validate_inputs(daily.drop(index=1))
        daily.loc[1, "date"] = daily.loc[0, "date"]
        with self.assertRaises(ValueError):
            depth.validate_inputs(daily)

    def test_unknown_station_and_incomplete_schema_are_rejected(self):
        daily = valid_input()
        with self.assertRaises(ValueError):
            depth.validate_inputs(daily.drop(columns="lo"))
        daily.loc[daily.station_id.eq(380), "station_id"] = 999
        with self.assertRaises(ValueError):
            depth.validate_inputs(daily)

    def test_nonfinite_curves_or_ranges_are_rejected(self):
        for field in ("SNOTEL", "iSnobal", "Noah-MP", "lo", "hi"):
            with self.subTest(field=field):
                daily = valid_input()
                daily.loc[0, field] = np.nan
                with self.assertRaises(ValueError):
                    depth.validate_inputs(daily)

    def test_range_order_and_center_membership_are_checked(self):
        for low, high, center in ((120., 100., 110.), (70., 110., 120.)):
            with self.subTest(values=(low, high, center)):
                daily = valid_input()
                daily.loc[0, ["lo", "hi", "iSnobal"]] = low, high, center
                with self.assertRaises(ValueError):
                    depth.validate_inputs(daily)

    def test_no_swe_timing_is_required_or_added(self):
        figure = depth.build_figure(valid_input())
        self.assertEqual(tuple(figure.get_size_inches()), (5.85, 3.4))
        self.assertEqual(len(figure.axes), 3)
        for ax in figure.axes:
            self.assertEqual(len(ax.lines), 3)
            self.assertEqual([line.get_marker() for line in ax.lines], ["None"] * 3)
            self.assertEqual(len(ax.collections), 1)
            self.assertEqual(len(ax.patches), 1)

    def test_all_daily_values_are_retained_despite_january_june_view(self):
        daily = valid_input().sample(frac=1, random_state=42)
        figure = depth.build_figure(daily)
        expected = pd.date_range("2025-10-01", "2026-06-30")
        for sid, ax in zip(depth.STATIONS, figure.axes):
            source = daily.loc[daily.station_id.eq(sid)].sort_values("date")
            for model, line in zip(depth.COLORS, ax.lines):
                np.testing.assert_array_equal(line.get_ydata(), source[model].to_numpy())
                self.assertEqual(len(line.get_xdata()), 273)
                self.assertTrue(pd.DatetimeIndex(line.get_xdata()).equals(expected))
            self.assertEqual(tuple(ax.get_xlim()), tuple(depth.mdates.date2num(pd.to_datetime(["2026-01-01", "2026-06-30"]))))


if __name__ == "__main__":
    unittest.main()
