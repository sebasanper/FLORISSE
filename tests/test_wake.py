# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import numpy as np
import unittest

import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData

import outputClasses.outputs

import florisCoreFunctions.windPlantFunctions as wPFs

class TestWake(unittest.TestCase):
    # TODO: Extend this with more checks
    def test_deflection_direction(self):
        # Generate and test a wake
        model = inputClasses.modelData.modelData(0, 0, 0)
        layout = layouts.Layout2(True)
        cSet = inputClasses.controlSettings.Neutral(layout)
        output = outputClasses.outputs.powerOutput(model, layout, cSet)

        # Set some required attributes in the output
        output.windSpeed[0] = 8
        output = wPFs.computeCpCtPoweraI(layout, cSet, output, 0)
        output.TI[0] = layout.TI_0
        wake = model.wake(model, layout, cSet, output, 0)

        self.assertTrue(hasattr(wake, 'V'))
        self.assertTrue(hasattr(wake, 'B'))
        self.assertTrue(hasattr(wake, 'displ'))

if __name__ == '__main__':
    unittest.main()
