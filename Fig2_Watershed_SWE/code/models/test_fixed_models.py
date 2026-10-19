"""Meaningful fixed-domain/weight/calendar regression tests for Fig.2."""
import unittest

import numpy as np
import pandas as pd

from aggregate_fig2_common_models import aggregate_fixed, check_calendar, coverage_precision_audit


class FixedModelsTest(unittest.TestCase):
    def setUp(self):
        self.w=np.array([[1.0,.25],[0,.75]])
        self.i=np.array([[[0,100],[np.nan,200]],[[200,0],[np.nan,0]]],dtype=float)
        self.n=np.zeros_like(self.i)

    def test_fractional_area_mean_including_finite_zero(self):
        i,n=aggregate_fixed(self.i,self.n,self.w)
        np.testing.assert_allclose(i,[87.5,100])
        np.testing.assert_allclose(n,[0,0])

    def test_equal_cell_weighting_would_change_result(self):
        i,_=aggregate_fixed(self.i,self.n,self.w)
        self.assertNotAlmostEqual(i[0],self.i[0][self.w>0].mean())

    def test_missing_retained_cell_aborts_instead_of_renormalizing(self):
        self.i[0,0,1]=np.nan
        with self.assertRaisesRegex(ValueError,"Nonfinite"):
            aggregate_fixed(self.i,self.n,self.w)

    def test_infinite_retained_cell_aborts(self):
        self.n[1,0,0]=np.inf
        with self.assertRaisesRegex(ValueError,"Nonfinite"):
            aggregate_fixed(self.i,self.n,self.w)

    def test_negative_swe_aborts(self):
        self.i[0,0,0]=-1
        with self.assertRaisesRegex(ValueError,"negative SWE"):
            aggregate_fixed(self.i,self.n,self.w)

    def test_invalid_outside_fixed_domain_ignored(self):
        self.n[:,1,0]=np.nan
        self.assertTrue(np.isfinite(aggregate_fixed(self.i,self.n,self.w)[0]).all())

    def test_invalid_weights_abort(self):
        self.w[0,1]=-1
        with self.assertRaisesRegex(ValueError,"nonnegative"):
            aggregate_fixed(self.i,self.n,self.w)

    def test_mask_cannot_add_cells_without_watershed_area(self):
        with self.assertRaisesRegex(ValueError,"fixed footprint"):
            aggregate_fixed(self.i,self.n,self.w,np.ones((2,2),bool))

    def test_shape_mismatch_aborts(self):
        with self.assertRaisesRegex(ValueError,"shape"):
            aggregate_fixed(self.i[:,0,:],self.n,self.w)

    def test_full_calendar_with_leap_days_passes(self):
        frames=[]
        for wy in range(2018,2027):
            end=f"{wy}-06-30" if wy==2026 else f"{wy}-09-30"
            frames.append(pd.DataFrame({"date":pd.date_range(f"{wy-1}-10-01",end),"wy":wy}))
        d=pd.concat(frames,ignore_index=True)
        self.assertEqual(len(check_calendar(d)),3195)
        without_leap=d[d.date.ne(pd.Timestamp("2020-02-29"))]
        with self.assertRaisesRegex(ValueError,"calendar"):
            check_calendar(without_leap)
        duplicated=pd.concat([d,d.tail(1)])
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            check_calendar(duplicated)

    def test_float32_coverage_omission_check(self):
        a=coverage_precision_audit(self.w.astype(np.float32))
        self.assertEqual(a["invisible_single_cell_omissions"],0)
        self.assertEqual(a["single_cell_omissions_tested"],3)


if __name__=="__main__":
    unittest.main()
