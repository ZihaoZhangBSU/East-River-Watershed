"""Synthetic operator tests with independently known overlap/weight results."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
from affine import Affine
from rasterio.crs import CRS
import rasterio

source=Path(__file__).resolve().parents[1]/'extract_fig2_common_aso.py'
spec=importlib.util.spec_from_file_location('aso',source)
aso=importlib.util.module_from_spec(spec);spec.loader.exec_module(aso)
CRS_M=CRS.from_epsg(32613)
SRC=Affine(1,0,500000,0,-1,4000000)
DST=Affine(2,0,500000,0,-2,4000000)

class RemapTests(unittest.TestCase):
    def test_missing_zero_extreme_and_fractional_weights(self):
        values=np.array([[0,2,np.nan,100000],[0,2,np.nan,100000]],dtype=float)
        swe,c,legacy,raw,pad=aso.remap(values,SRC,CRS_M,DST,CRS_M,(1,2))
        np.testing.assert_allclose(swe,[[1,100000]],rtol=0,atol=1e-6)
        np.testing.assert_allclose(c,[[1,.5]],rtol=0,atol=1e-6)
        mean,fraction,mask=aso.basin_mean(swe,np.array([[1,.25]]),np.ones((1,2),bool))
        self.assertAlmostEqual(mean,20000.8)
        self.assertEqual(fraction,1.)
        # Multiplying by C_j would give a different mean; it is not the operator.
        self.assertGreater(abs(mean-(1+12500)/1.125),8000)

    def test_changing_missing_mask_has_explicit_coverage(self):
        values=np.array([[0,2,np.nan,np.nan],[0,2,np.nan,np.nan]])
        swe,c,*_=aso.remap(values,SRC,CRS_M,DST,CRS_M,(1,2))
        mean,fraction,mask=aso.basin_mean(swe,np.array([[1,.25]]),np.ones((1,2),bool))
        self.assertEqual(mean,1.)
        self.assertEqual(fraction,.8)
        np.testing.assert_array_equal(mask,[[True,False]])

    def test_source_extent_requires_full_target_denominator(self):
        target=Affine(2,0,499999.25,0,-2,4000000)
        swe,c,legacy,raw,pad=aso.remap(np.ones((2,2))*12,SRC,CRS_M,target,CRS_M,(1,2))
        np.testing.assert_allclose(legacy,[[1,1]],rtol=0,atol=1e-6)
        np.testing.assert_allclose(c,[[.625,.375]],rtol=0,atol=1e-6)
        self.assertEqual(swe[0,0],12)
        self.assertTrue(np.isnan(swe[0,1]))

    def test_threshold_includes_exact_half(self):
        target=Affine(2,0,499999,0,-2,4000000)
        swe,c,*_=aso.remap(np.ones((2,2))*7,SRC,CRS_M,target,CRS_M,(1,2))
        np.testing.assert_allclose(c,[[.5,.5]],rtol=0,atol=1e-6)
        np.testing.assert_allclose(swe,[[7,7]])

    def test_fixed_domain_does_not_expand_to_observations(self):
        mean,fraction,valid=aso.basin_mean(np.array([[0.,900.]]),np.array([[.25,1.]]),np.array([[True,False]]))
        self.assertEqual(mean,0.)
        self.assertEqual(fraction,1.)
        np.testing.assert_array_equal(valid,[[True,False]])

    def test_negative_values_are_not_screened(self):
        swe,c,*_=aso.remap(np.array([[-1000,0],[-1000,0]]),SRC,CRS_M,DST,CRS_M,(1,1))
        self.assertEqual(swe[0,0],-500.)
        self.assertEqual(c[0,0],1.)

    def test_fractional_boundary_area_not_equal_cell_mean(self):
        mean,fraction,_=aso.basin_mean(np.array([[0.,100.]]),np.array([[.25,1.]]),np.ones((1,2),bool))
        self.assertEqual(mean,80.)
        self.assertEqual(fraction,1.)

    def test_north_south_target_orientation(self):
        northup,cn,*_=aso.remap(np.arange(16).reshape(4,4),SRC,CRS_M,DST,CRS_M,(2,2))
        southup,cs,*_=aso.remap(np.arange(16).reshape(4,4),SRC,CRS_M,Affine(2,0,500000,0,2,3999996),CRS_M,(2,2))
        np.testing.assert_allclose(northup,southup[::-1])
        np.testing.assert_allclose(cn,cs[::-1])

    def test_source_mask_units_scale_offset_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            source=root/'input.tif'
            # Declared -9999 NoData must never enter averaging. Valid zero and
            # high values are transformed by GDAL scale/offset before averaging.
            with rasterio.open(source,'w',driver='GTiff',height=2,width=4,count=1,dtype='float32',
                               crs=CRS_M,transform=SRC,nodata=-9999) as dst:
                dst.write(np.array([[0,2,-9999,100000],[0,2,-9999,100000]],'float32'),1)
                dst.scales=(2.,);dst.offsets=(3.,)
                dst.update_tags(description='Snow water equivalent [mm]')
            flights=root/'flights.csv';flights.write_text('date,source_file\n2020-04-01,input.tif\n')
            grid=root/'grid.npz'
            np.savez(grid,weights=np.array([[1,.25]]),fixed_domain_mask=np.ones((1,2),bool),
                     transform=np.array(tuple(DST)),crs=CRS_M.to_wkt())
            rows,manifest=aso.extract(root,flights,grid,root/'out')
            self.assertAlmostEqual(rows[0]['aso_watershed_mean_mm'],40004.6)
            self.assertEqual(rows[0]['aso_valid_area_fraction'],1.)
            self.assertEqual(manifest['sources'][0]['scale'],2.)
            self.assertEqual(manifest['sources'][0]['offset'],3.)
            self.assertEqual(aso.sha256(source),manifest['sources'][0]['sha256'])

    def test_grid_must_declare_domain_and_geometry(self):
        with tempfile.TemporaryDirectory() as directory:
            grid=Path(directory)/'bad.npz'
            np.savez(grid,weights=np.ones((1,2)))
            with self.assertRaisesRegex(ValueError,'must include'):
                aso.load_grid(grid)

if __name__=='__main__':
    unittest.main()
