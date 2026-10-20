"""Synthetic checks; no research observations are distributed with these tests."""
from pathlib import Path
import sys
import tempfile
import unittest

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Transformer
from shapely.geometry import Polygon

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import Supp_ASO_Extremes as plot


def fixture():
    dates = [d for d, n in plot.THRESHOLD_COUNTS.items() for _ in range(n)]
    x = 320000. + 50 * np.arange(78)
    y = np.full(78, 4310000.)
    lon, lat = Transformer.from_crs(32613, 4326, always_xy=True).transform(x, y)
    return pd.DataFrame({"date": dates, "row": range(78), "column": range(78),
        "x_m": x, "y_m": y, "longitude": lon, "latitude": lat,
        "crs": "EPSG:32613", "swe_m": 5.1, "snow_depth_m": 12.})


class ProximityTests(unittest.TestCase):
    def test_cross_flight_inclusive_125_m(self):
        q = pd.DataFrame({"x_m": [0., 125., 250.0001, 0.], "y_m": [0.] * 4,
                          "date": ["a", "b", "c", "a"]})
        np.testing.assert_array_equal(plot.recurrence_mask(q), [True, True, False, True])

    def test_same_flight_is_never_recurrence(self):
        q = pd.DataFrame({"x_m": [0., 0.], "y_m": [0., 0.], "date": ["a", "a"]})
        self.assertFalse(plot.recurrence_mask(q).any())


class BoundaryTests(unittest.TestCase):
    def test_real_crs_transformation_and_feature_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.geojson"
            ring = [(-107., 38.7), (-106.9, 38.7), (-106.9, 38.8), (-107., 38.7)]
            gpd.GeoDataFrame({"huc10": ["1402000102", "0000000000"]},
                geometry=[Polygon(ring), Polygon(ring)], crs=4269).to_file(path, driver="GeoJSON")
            boundary, metadata = plot.load_boundary(path)
            x, y = Transformer.from_crs(4269, 32613, always_xy=True).transform(*ring[0])
            self.assertAlmostEqual(boundary.exterior.coords[0][0], x, places=6)
            self.assertAlmostEqual(boundary.exterior.coords[0][1], y, places=6)
            self.assertEqual(metadata["huc10"], "1402000102")
            self.assertEqual(metadata["source_vertex_count"], 4)

    def test_missing_watershed_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.geojson"
            gpd.GeoDataFrame({"huc10": ["0000000000"]},
                geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 0)])],
                crs=4326).to_file(path, driver="GeoJSON")
            with self.assertRaisesRegex(ValueError, "exactly one"):
                plot.load_boundary(path)


class InputAndRenderingTests(unittest.TestCase):
    def test_wrong_zero_count_flight_and_wrong_order_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "flights.csv"
            dates = list(plot.REFERENCE_FLIGHTS)
            pd.DataFrame({"date": dates}).to_csv(path, index=False)
            self.assertEqual(plot.load_flights(path), dates)
            wrong = dates.copy()
            wrong[0] = "2018-04-01"  # Distinct date; original zero-count flight lost.
            pd.DataFrame({"date": wrong}).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "exact 14"):
                plot.load_flights(path)
            pd.DataFrame({"date": dates[::-1]}).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "chronological"):
                plot.load_flights(path)

    def test_no_source_mutation_or_lost_samples(self):
        q = fixture()
        prior = q.copy(deep=True)
        flights = list(plot.THRESHOLD_COUNTS) + [f"2025-01-{i:02d}" for i in range(1, 9)]
        boundary = Polygon([(315000, 4290000), (350000, 4290000),
                            (350000, 4330000), (315000, 4290000)])
        figure = plot.make_figure(q, flights, boundary)
        collections = figure.axes[0].collections
        self.assertEqual(sum(len(c.get_offsets()) for c in collections[:6]), 78)
        np.testing.assert_array_equal([r.get_width() for r in figure.axes[1].patches],
                                      list(plot.THRESHOLD_COUNTS.values()) + [0] * 8)
        pd.testing.assert_frame_equal(q, prior)
        plt.close(figure)

    def test_invalid_reference_records_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.csv"
            valid = fixture()
            valid.to_csv(path, index=False)
            self.assertEqual(len(plot.load_locations(path)), 78)
            bad = valid.copy()
            bad.loc[0, "swe_m"] = 5.0
            bad.loc[0, "snow_depth_m"] = 15.
            bad.to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "criterion"):
                plot.load_locations(path)
            bad = valid.copy()
            bad.loc[0, "longitude"] += .01
            bad.to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "coordinates disagree"):
                plot.load_locations(path)


if __name__ == "__main__":
    unittest.main()
