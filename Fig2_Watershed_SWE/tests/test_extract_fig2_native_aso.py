"""Scientific regression checks for native footprint means and coverage."""
import sys
sys.dont_write_bytecode = True

import unittest
from pathlib import Path
import numpy as np
from pyproj import CRS
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from shapely.geometry import box

from extract_fig2_native_aso import (calculate_mean, overlap_weights,
                                     projected_area_factor, swe_units_factor, validate_output_paths)


class NativeFootprintTests(unittest.TestCase):
    def test_fractional_boundary_nodata_nonfinite_unfiltered_and_full_denominator(self):
        # The polygon extends 1 m beyond the source raster, and bisects its first
        # column. Retain zero, a negative value, and a very large positive value.
        data = np.array([[1., 0., 4000.], [-0.5, -9999., np.nan]], dtype="float32")
        with MemoryFile() as mem:
            with mem.open(driver="GTiff", width=3, height=2, count=1, dtype="float32",
                          transform=from_origin(0, 2, 1, 1), crs="EPSG:32613", nodata=-9999.) as src:
                src.write(data, 1)
                src.update_tags(Band_1="ASO Snow Water Equivalent [m]")
                window, areas, watershed_area = overlap_weights(src, box(.5, 0, 4, 2), CRS("EPSG:32613"))
                factor, _ = swe_units_factor(src, "auto")
                values = src.read(1, window=window, masked=True).astype(float).filled(np.nan) * factor
                mean, fraction, valid_area = calculate_mean(values, areas, watershed_area)
        np.testing.assert_allclose(areas, [[.5, 1, 1], [.5, 1, 1]])
        self.assertEqual(watershed_area, 7.)
        self.assertEqual(valid_area, 3.)
        self.assertAlmostEqual(fraction, 3. / 7.)
        self.assertAlmostEqual(mean, (500. + 4000000. - 250.) / 3.)

    def test_geographic_degree_area_rejected(self):
        with self.assertRaisesRegex(ValueError, "geographic degree areas"):
            projected_area_factor(CRS("EPSG:4326"))

    def test_projected_feet_area_is_converted_to_square_metres(self):
        crs = CRS("EPSG:2227")
        self.assertAlmostEqual(projected_area_factor(crs), (1200. / 3937.) ** 2)

    def test_no_valid_pixels_fails(self):
        with self.assertRaisesRegex(ValueError, "No finite native"):
            calculate_mean(np.array([[np.nan]]), np.array([[1.]]), 1.)

    def test_unknown_units_fail(self):
        with MemoryFile() as mem:
            with mem.open(driver="GTiff", width=1, height=1, count=1, dtype="float32",
                          transform=from_origin(0, 2, 1, 1), crs="EPSG:32613") as src:
                with self.assertRaisesRegex(ValueError, "Cannot establish SWE units"):
                    swe_units_factor(src, "auto")
                self.assertEqual(swe_units_factor(src, "mm")[0], 1.)
                src.update_tags(pixel_size="50 m")
                with self.assertRaisesRegex(ValueError, "Cannot establish SWE units"):
                    swe_units_factor(src, "auto")

    def test_conflicting_swe_units_fail(self):
        with MemoryFile() as mem:
            with mem.open(driver="GTiff", width=1, height=1, count=1, dtype="float32",
                          transform=from_origin(0, 2, 1, 1), crs="EPSG:32613") as src:
                src.set_band_unit(1, "m")
                src.update_tags(Band_1="ASO Snow Water Equivalent [mm]")
                with self.assertRaisesRegex(ValueError, "Conflicting SWE unit"):
                    swe_units_factor(src, "auto")

    def test_output_paths_protect_inputs(self):
        with self.assertRaisesRegex(ValueError, "must be distinct"):
            validate_output_paths(Path("output.csv"), Path("./output.csv"), [])
        with self.assertRaisesRegex(ValueError, "must not overwrite"):
            validate_output_paths(Path("./source.tif"), Path("manifest.json"), [Path("source.tif")])
        with self.assertRaisesRegex(ValueError, "must not overwrite"):
            validate_output_paths(Path("output.csv"), Path("./flights.csv"), [Path("flights.csv")])


if __name__ == "__main__":
    unittest.main(verbosity=2)
