# -*- coding: utf-8 -*-
import imp
import unittest
import pickle

from Turbines.NREL5MW.NREL5MW import NREL5MWTurbine


class TestNREL5MWProperties(unittest.TestCase):
    """Tests for `FLORISSE\florisCoreFunctions\windPlant.py`."""
    def test_Attributes_of_NREL5MW_with_Pitch(self):
        # Select wind farm layout and specify controlset
        turbine = NREL5MWTurbine(True)

        self.assertTrue(turbine.usePitch)
        self.assertTrue(hasattr(turbine, 'betaLims'))

        self.assertEqual(turbine.rotorDiameter, 126)
        self.assertEqual(turbine.hubHeight, 90)
        self.assertEqual(turbine.NumBlades, 3)
        self.assertEqual(turbine.pP, 1.88)
        self.assertEqual(turbine.pT, 2.07)
        self.assertEqual(turbine.gE, 1.0)
        self.assertEqual(turbine.eta, 0.768)

        CpCtWithPitch = pickle.load(open('testData/CpCtWithPitch.p', 'rb'))
        for pitch, WSdict in CpCtWithPitch.items():
            for windSpeed, value in WSdict.items():
                self.assertAlmostEqual(turbine.Ct(windSpeed, pitch), value, 8)

    def test_Attributes_of_NREL5MW_without_Pitch(self):
        # Select wind farm layout and specify controlset
        turbine = NREL5MWTurbine(False)

        self.assertFalse(turbine.usePitch)
        self.assertFalse(hasattr(turbine, 'betaLims'))

        self.assertEqual(turbine.rotorDiameter, 126)
        self.assertEqual(turbine.hubHeight, 90)
        self.assertEqual(turbine.NumBlades, 3)
        self.assertEqual(turbine.pP, 1.88)
        self.assertEqual(turbine.pT, 2.07)
        self.assertEqual(turbine.gE, 1.0)
        self.assertEqual(turbine.eta, 0.768)

        CpCtWithoutPitch = pickle.load(open('testData/CpCtWithoutPitch.p', 'rb'))
        for windSpeed, value in CpCtWithoutPitch.items():
            self.assertAlmostEqual(turbine.Ct(windSpeed), value, 8)


if __name__ == '__main__':
    unittest.main() 
