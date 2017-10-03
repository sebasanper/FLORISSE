# -*- coding: utf-8 -*-
import numpy as np
import unittest

import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData

import outputClasses.outputs

import florisCoreFunctions.windPlantFunctions as wPFs
from florisCoreFunctions.wake import Wake


class TestWake(unittest.TestCase):
    # TODO, extend this with more cases
    def setUp(self):
        # Generate and test a wake
        model = inputClasses.modelData.modelData(0, 0, 0)
        layout = layouts.Layout2(True)
        cSet = inputClasses.controlSettings.Neutral(layout)
        output = outputClasses.outputs.powerOutput(model, layout, cSet)

        # Set some required attributes in the output
        output.windSpeed[0] = 8
        output = wPFs.computeCpCtPoweraI(layout, cSet, output, 0)
        output.TI[0] = layout.TI_0
        self.wake = Wake(model, layout, cSet, output, 0)

    def test_deflection_direction(self):

        self.assertTrue(hasattr(self.wake, 'displ'))
        sh = np.zeros([1, 1])
        np.testing.assert_array_almost_equal(self.wake.V(8, -10, sh, sh), 8, 5)

if __name__ == '__main__':
    unittest.main()
