"""Depth-specific masks, weighting, complete calendars and temporal statistics."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

code=Path(__file__).resolve().parents[1]
def module(name):
    spec=importlib.util.spec_from_file_location(name,code/(name+'.py'))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result
extract=module('extract_watershed_depth');plotter=module('Supp_Watershed_Depth')

class DepthOperatorTests(unittest.TestCase):
    def test_fractional_watershed_weights_are_not_equal_cell_weights(self):
        mean,fraction,mask,error=extract.weighted_mean(np.array([[0.,100.]]),np.array([[.25,.75]]),np.ones((1,2),bool))
        self.assertEqual(mean,75.);self.assertEqual(fraction,1.)

    def test_finite_zero_negative_and_extreme_values_remain(self):
        values=np.array([[0.,-10.,20000.]])
        mean,*_=extract.weighted_mean(values,np.ones((1,3)),np.ones((1,3),bool))
        self.assertEqual(mean,(20000.-10.)/3)

    def test_missing_daily_depth_aborts_instead_of_renormalizing(self):
        with self.assertRaisesRegex(ValueError,'Incomplete finite model DEPTH'):
            extract.weighted_mean(np.array([[0.,np.nan]]),np.ones((1,2)),np.ones((1,2),bool),True)

    def test_ASO_uses_own_depth_finite_support_and_keeps_zero(self):
        mean,fraction,mask,error=extract.depth_from_corrected(np.array([[0.,100.,np.nan]]),
            np.array([[.8,.6,.4]]),np.array([[.25,.25,.5]]),np.ones((1,3),bool))
        self.assertEqual(mean,50.);self.assertEqual(fraction,.5)
        np.testing.assert_array_equal(mask,np.array([[True,True,False]]))

    def test_half_coverage_is_inclusive_and_no_Cj_multiplier_is_used(self):
        mean,*_=extract.depth_from_corrected(np.array([[0.,100.]]),np.array([[1.,.5]]),
            np.ones((1,2)),np.ones((1,2),bool))
        self.assertEqual(mean,50.)

    def test_finite_depth_below_own_coverage_threshold_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'Cj>=0.5'):
            extract.depth_from_corrected(np.array([[0.]]),np.array([[.499]]),np.ones((1,1)),np.ones((1,1),bool))

    def test_model_and_ASO_do_not_silently_drop_positive_watershed_targets(self):
        with self.assertRaises(ValueError):
            extract.weighted_mean(np.ones((1,2)),np.ones((1,2)),np.array([[True,False]]))

    def test_leap_days_and_partial_2026_calendar(self):
        self.assertEqual(len(extract.calendar(2020)),366)
        self.assertIn(pd.Timestamp('2020-02-29'),extract.calendar(2020))
        self.assertEqual(len(extract.calendar(2026)),273)
        self.assertEqual(str(extract.calendar(2026)[-1].date()),'2026-06-30')

class DepthInputAndStatisticsTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        (self.root/'models').mkdir();(self.root/'aso').mkdir()
        dates=pd.DatetimeIndex([])
        for wy in plotter.YEARS:dates=dates.append(extract.calendar(wy))
        self.model=pd.DataFrame({'date':dates,'wy':dates.year+(dates.month>=10),
            'variable':'Depth','depth_i_mm':np.arange(len(dates),dtype=float),
            'depth_n_mm':np.arange(len(dates),dtype=float)+1,
            'fixed_model_area_fraction':1.,'daily_common_valid_area_fraction':1.,
            'fixed_model_cells':842,'state_hour_utc':23})
        aso_dates=pd.to_datetime(plotter.FLIGHT_DATES);area=748679519.6193271
        self.aso=pd.DataFrame({'date':aso_dates,'wy':aso_dates.year,'variable':'Depth',
            'aso_watershed_mean_mm':100.,'aso_valid_area_fraction':.9,
            'retained_target_watershed_area_m2':.9*area,'fixed_model_watershed_area_m2':area})

    def tearDown(self):plotter.plt.close('all');self.temp.cleanup()
    def load(self):
        self.model.to_csv(self.root/'models/watershed_depth_daily.csv',index=False)
        self.aso.to_csv(self.root/'aso/aso_depth_flights.csv',index=False)
        return plotter.load_inputs(self.root)

    def test_complete_depth_calendar_and_flights_accepted(self):
        model,aso=self.load();self.assertEqual(len(model),3195);self.assertEqual(len(aso),14)

    def test_SWE_label_cannot_be_replayed_as_depth(self):
        self.aso.loc[0,'variable']='SWE'
        with self.assertRaisesRegex(ValueError,'identify Depth'):self.load()

    def test_calendar_gap_is_rejected(self):
        self.model=self.model.drop(index=100)
        with self.assertRaisesRegex(ValueError,'every displayed day'):self.load()

    def test_missing_ASO_flight_is_rejected(self):
        self.aso=self.aso.iloc[:-1]
        with self.assertRaisesRegex(ValueError,'14 documented'):self.load()

    def test_nonfinite_daily_depth_is_rejected(self):
        self.model.loc[0,'depth_i_mm']=np.nan
        with self.assertRaisesRegex(ValueError,'must be finite'):self.load()

    def test_partial_daily_model_support_is_rejected(self):
        self.model.loc[0,'daily_common_valid_area_fraction']=.99
        with self.assertRaisesRegex(ValueError,'complete finite support'):self.load()

    def test_inconsistent_ASO_depth_area_denominator_is_rejected(self):
        self.aso.loc[0,'aso_valid_area_fraction']=.8
        with self.assertRaisesRegex(ValueError,'retained / fixed'):self.load()

    def test_full_historical_periods_and_2026_partial_are_used(self):
        stats=plotter.model_statistics(*self.load())
        self.assertEqual(stats.n_dates.tolist(),[365,365,366,365,365,365,366,365,273])
        self.assertEqual(stats.end_date.iloc[0],'2018-09-30')
        self.assertEqual(stats.end_date.iloc[-1],'2026-06-30')

    def test_temporal_r_bias_sign_and_relative_denominator(self):
        stats=plotter.paired_model_statistics(np.array([0.,1.,2.]),np.array([1.,2.,3.]))
        self.assertAlmostEqual(stats['pearson_r'],1.)
        self.assertEqual(stats['mean_bias_noah_minus_isnobal_mm'],1.)
        self.assertEqual(stats['relative_difference_noah_vs_isnobal_pct'],100.)

    def test_zero_denominator_and_constant_correlation_are_NA_not_zero(self):
        stats=plotter.paired_model_statistics(np.zeros(3),np.ones(3))
        self.assertTrue(np.isnan(stats['pearson_r']))
        self.assertTrue(np.isnan(stats['relative_difference_noah_vs_isnobal_pct']))
        self.assertIn('r = NA',plotter.statistics_annotation(pd.Series(stats)))
        self.assertIn('rMB = NA',plotter.statistics_annotation(pd.Series(stats)))

    def test_incomplete_temporal_pairs_do_not_change_support_silently(self):
        with self.assertRaises(ValueError):plotter.paired_model_statistics(np.array([0.,np.nan]),np.ones(2))

    def test_nine_panels_eighteen_curves_and_fourteen_depth_markers(self):
        model,aso=self.load();figures=[]
        with patch.object(plotter.plt.Figure,'savefig',lambda fig,*a,**k:figures.append(fig)):
            paths=plotter.render(model,aso,self.root,None)
        fig=figures[0];self.assertEqual(len(fig.axes),9)
        self.assertEqual(sum(len(ax.lines) for ax in fig.axes),18)
        self.assertEqual(sum(len(c.get_offsets()) for ax in fig.axes for c in ax.collections),14)
        self.assertEqual([p.name for p in paths],['JHM_S_watershed_depth.pdf','JHM_S_watershed_depth.png'])
        self.assertEqual(fig._supylabel.get_text(),'Watershed snow depth (mm)')
        self.assertEqual(tuple(fig.get_size_inches()),(6.5,5.5))
        self.assertEqual(fig.axes[-1].get_title(),'WY2026 (Oct-Jun)')

if __name__=='__main__':unittest.main()
