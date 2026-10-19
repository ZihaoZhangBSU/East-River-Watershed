"""Exercise scientifically consequential corrupted-input cases for Fig.2.

Fixtures are synthetic calendar/coverage records only; no research data is
embedded. Zero-SWE fixtures ensure valid snow-free states stay admissible.
"""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

_path=Path(__file__).resolve().parents[1]/"Fig2_Watershed_SWE.py"
_spec=importlib.util.spec_from_file_location("fig2_generator",_path)
fig2=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig2)


class FigureInputValidationTests(unittest.TestCase):
    def setUp(self):
        self.scratch=Path.cwd()/"work"/"fig2_input_tests"
        self.scratch.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(prefix="fig2_input_test_",dir=self.scratch)
        self.root=Path(self.temp.name)
        (self.root/"models").mkdir()
        (self.root/"aso").mkdir()
        dates=pd.date_range("2017-10-01","2026-06-30")
        self.models=pd.DataFrame({"date":dates,
            "wy":dates.year+(dates.month>=10).astype(int),
            "swe_i_mm":np.zeros(len(dates)),"swe_n_mm":np.zeros(len(dates)),
            "fixed_model_area_fraction":1.,"daily_common_valid_area_fraction":1.,
            "fixed_model_cells":842,"state_hour_utc":23})
        flights=pd.to_datetime(fig2.FLIGHT_DATES)
        basin=748679519.6193271
        self.aso=pd.DataFrame({"date":flights,
            "wy":flights.year+(flights.month>=10).astype(int),
            "aso_watershed_mean_mm":np.zeros(len(flights)),
            "aso_valid_area_fraction":.75,
            "retained_target_watershed_area_m2":.75*basin,
            "fixed_model_watershed_area_m2":basin})

    def tearDown(self):
        if not self.root.resolve().is_relative_to(self.scratch.resolve()):
            raise AssertionError("Refusing cleanup outside the test scratch directory")
        self.temp.cleanup()

    def load(self):
        self.models.to_csv(self.root/"models/Fig2_common_grid_daily.csv",index=False)
        self.aso.to_csv(self.root/"aso/aso_common_grid_flights.csv",index=False)
        return fig2.load_inputs(self.root)

    def test_complete_calendar_and_finite_zeros_are_accepted(self):
        models,aso=self.load()
        self.assertEqual(len(models),3195)
        self.assertEqual(len(aso),14)
        self.assertEqual(float(models.swe_i_mm.sum()),0.)
        self.assertEqual(float(aso.aso_watershed_mean_mm.sum()),0.)

    def test_unfiltered_supplied_finite_aso_values_are_not_clipped(self):
        self.aso.loc[0,"aso_watershed_mean_mm"]=-.1
        _,aso=self.load()
        self.assertEqual(float(aso.aso_watershed_mean_mm.iloc[0]),-.1)

    def test_missing_ordinary_day_rejected(self):
        self.models=self.models[self.models.date.ne(pd.Timestamp("2021-07-11"))]
        with self.assertRaisesRegex(ValueError,"every displayed day"):
            self.load()

    def test_missing_leap_day_rejected(self):
        self.models=self.models[self.models.date.ne(pd.Timestamp("2020-02-29"))]
        with self.assertRaisesRegex(ValueError,"leap days"):
            self.load()

    def test_duplicate_flight_rejected(self):
        self.aso=pd.concat([self.aso,self.aso.iloc[[0]]],ignore_index=True)
        with self.assertRaisesRegex(ValueError,"unique"):
            self.load()

    def test_wrong_model_water_year_rejected(self):
        self.models.loc[0,"wy"]=2017
        with self.assertRaisesRegex(ValueError,"water-year labels"):
            self.load()

    def test_wrong_flight_water_year_rejected(self):
        self.aso.loc[0,"wy"]=2019
        with self.assertRaisesRegex(ValueError,"water-year labels"):
            self.load()

    def test_missing_model_support_rejected_even_with_finite_swe(self):
        self.models.loc[0,"daily_common_valid_area_fraction"]=.999
        with self.assertRaisesRegex(ValueError,"complete finite support"):
            self.load()

    def test_mislabeled_flight_coverage_rejected(self):
        self.aso.loc[0,"aso_valid_area_fraction"]=.9
        with self.assertRaisesRegex(ValueError,"retained / fixed"):
            self.load()

    def test_nonfinite_swe_rejected(self):
        self.models.loc[0,"swe_i_mm"]=np.nan
        with self.assertRaisesRegex(ValueError,"finite and nonnegative"):
            self.load()

    def test_changed_flight_coverage_denominator_rejected(self):
        # Keep the ratio unchanged to expose an inconsistent fixed-area basis.
        self.aso.loc[0,"fixed_model_watershed_area_m2"]*=2
        self.aso.loc[0,"retained_target_watershed_area_m2"]*=2
        with self.assertRaisesRegex(ValueError,"fixed.*area|area.*fixed"):
            self.load()


if __name__=="__main__":
    unittest.main()
