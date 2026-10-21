"""Meaningful synthetic tests for equal-cell model-specific F5 scores."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
source=Path(__file__).resolve().parents[1]/'Fig5_ASO_Performance.py'
spec=importlib.util.spec_from_file_location('f5',source)
f5=importlib.util.module_from_spec(spec);spec.loader.exec_module(f5)

class ScoreTests(unittest.TestCase):
    def test_model_specific_pairs(self):
        observation=np.array([[0.,2.,np.nan,100.]])
        domain=np.ones(observation.shape,bool)
        i=f5.scores(observation,np.array([[1.,2.,100.,np.nan]]),domain)
        n=f5.scores(observation,np.array([[0.,np.nan,0.,100.]]),domain)
        self.assertEqual(i['n_cells'],2)
        self.assertEqual(n['n_cells'],2)
        self.assertEqual(i['aso_mean_mm'],1.)
        self.assertEqual(n['aso_mean_mm'],50.)
        self.assertEqual(i['mb_mm'],.5)
        self.assertAlmostEqual(i['rmse_mm'],np.sqrt(.5))
        self.assertEqual(n['rmse_mm'],0.)

    def test_equal_cell_score_preserves_zero(self):
        # Basin fractional weights may differ; F5 gives each eligible target
        # one vote. It must not substitute an area-weighted basin mean.
        observation=np.array([[0.,100.]])
        result=f5.scores(observation,np.array([[100.,100.]]),np.ones((1,2),bool))
        self.assertEqual(result['n_cells'],2)
        self.assertEqual(result['aso_mean_mm'],50.)
        self.assertEqual(result['mb_mm'],50.)
        self.assertAlmostEqual(result['rmse_mm'],np.sqrt(5000.))

    def test_domain_and_missing_values_are_distinct(self):
        observation=np.array([[0.,100.,1000.]])
        model=np.array([[0.,np.nan,-1000.]])
        result=f5.scores(observation,model,np.array([[True,True,False]]))
        self.assertEqual(result['n_cells'],1)
        self.assertEqual(result['rmse_mm'],0.)

    def test_elevation_boundaries(self):
        elevation=np.array([2895.9,2896.,3352.9,3353.])
        counts=[]
        for name in ['Lower','Middle','Upper']:
            lo,hi=f5.BANDS[name]
            mask=np.ones(elevation.shape,bool)
            if lo is not None:mask &= elevation>=lo
            if hi is not None:mask &= elevation<hi
            counts.append(int(mask.sum()))
        self.assertEqual(counts,[1,2,1])

    def test_spatial_pearson_retains_sign_instead_of_r_squared(self):
        observed=np.array([[0.,1.,2.]])
        result=f5.scores(observed,np.array([[2.,1.,0.]]),np.ones((1,3),bool))
        self.assertAlmostEqual(result['r'],-1.)

    def test_constant_field_correlation_is_undefined_not_zero(self):
        domain=np.ones((1,3),bool)
        self.assertIsNone(f5.scores(np.array([[0.,1.,2.]]),np.ones((1,3)),domain)['r'])
        self.assertIsNone(f5.scores(np.ones((1,3)),np.array([[0.,1.,2.]]),domain)['r'])

    def test_nonconstant_zero_spatial_correlation_is_valid(self):
        result=f5.scores(np.array([[-1.,0.,1.]]),np.array([[1.,-2.,1.]]),np.ones((1,3),bool))
        self.assertEqual(result['r'],0.)

if __name__=='__main__':
    unittest.main()
