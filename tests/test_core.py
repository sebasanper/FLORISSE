# -*- coding: utf-8 -*-
import imp
import unittest
import pickle

import florisCoreFunctions.windPlant as windPlant
import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData


class TestFlorisCoreWithGaussian(unittest.TestCase):

    def setUp(self):
        # Select a velocity, deflection and wake summing model
        self.model = inputClasses.modelData.modelData(2, 1, 2)

    """Tests for `FLORISSE\florisCoreFunctions\windPlant.py`."""
    def test_power_of_single_turbine_layout_with_NREL5MW_Pitch(self):
        # Select wind farm layout and specify controlset
        layout = layouts.Layout1(True)
        cSet = inputClasses.controlSettings.Neutral(layout)

        # Run the model and get an output object
        outputNeutral = windPlant.windPlant(self.model, layout, cSet, False)

        """Does a full model run give the correct anwer?"""
        self.assertAlmostEqual(sum(outputNeutral.power)[0], 1162592.50473, 4)

    def test_power_of_2_by_2_layout_with_NREL5MW_Pitch(self):
        # Select wind farm layout and specify controlset
        layout = layouts.Layout2(True)
        cSet = inputClasses.controlSettings.Neutral(layout)

        # In case of model changes, create new data with
        # powerDataWindRose2by2[deg] = sum(outputData['powerOut'])
        # pickle.dump(powerDataWindRose2by2, open('powerDataWindRose2by2.p', 'wb'))
        powerDataWindRose2by2 = pickle.load(open('testData/powerDataWindRose2by2.p', 'rb'))
        for key, value in powerDataWindRose2by2.items():
            layout.windDirection = key
            # Run the model and get a new output object
            outputNeutral = windPlant.windPlant(self.model, layout, cSet, False)

            """Does the model still behave as expected"""
            self.assertAlmostEqual(sum(outputNeutral.power)[0], value, 8)

if __name__ == '__main__':
    unittest.main()
