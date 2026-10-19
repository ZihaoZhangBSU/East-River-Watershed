"""Scientific support-policy checks for the narrow S2 coverage correction."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

_path=Path(__file__).resolve().parents[1]/"code/S2_Three_Source_Density.py"
_spec=importlib.util.spec_from_file_location("s2_density",_path)
s2=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s2)


class DensitySupportTests(unittest.TestCase):
    def setUp(self):
        self.fields={"basin":np.ones((1,3),bool)}
        for source in ("ASO","iSnobal","Noah-MP"):
            self.fields[source+"_SWE"]=np.array([[0.,30.,50.]])
            self.fields[source+"_depth"]=np.array([[200.,200.,200.]])

    def test_finite_zero_swe_is_retained(self):
        self.assertTrue(s2.support(self.fields)[0,0])

    def test_strict_100mm_aso_depth_screen(self):
        self.fields["ASO_depth"][0]=[99.9,100.,100.1]
        np.testing.assert_array_equal(s2.support(self.fields),[[False,False,True]])

    def test_model_missed_snow_excluded_from_all_three(self):
        self.fields["iSnobal_depth"][0,0]=0
        self.fields["Noah-MP_depth"][0,1]=0
        np.testing.assert_array_equal(s2.support(self.fields),[[False,False,True]])

    def test_any_missing_component_excludes_shared_sample(self):
        self.fields["ASO_SWE"][0,0]=np.nan
        self.fields["Noah-MP_depth"][0,1]=np.nan
        np.testing.assert_array_equal(s2.support(self.fields),[[False,False,True]])

    def test_supplied_negative_swe_is_not_density_clipped(self):
        self.fields["ASO_SWE"][0,0]=-.1
        self.assertTrue(s2.support(self.fields)[0,0])

    def test_cells_outside_watershed_support_excluded(self):
        self.fields["basin"][0,0]=False
        np.testing.assert_array_equal(s2.support(self.fields),[[False,True,True]])


if __name__=="__main__":
    unittest.main()
