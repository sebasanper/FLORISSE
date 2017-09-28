# -*- coding: utf-8 -*-
import numpy as np
import unittest

import inputClasses.layouts as layouts
import inputClasses.controlSettings


class TestControlSet(unittest.TestCase):

    """Tests for `FLORISSE\florisCoreFunctions\windPlant.py`."""
    def test_neutral_controlSet(self):
        layout = layouts.Layout1(True)
        cSet = inputClasses.controlSettings.Neutral(layout)
        for angle in np.arange(30, -30, -5):
            # TODO: don't use abs() but sign derived from wakeDir
            cSet.yawAngles[0] = angle
            self.assertAlmostEqual(abs(np.radians(angle)), cSet.alphas[0], 8,
                                   'yaw test filled at angle = %f'%angle)
            cSet.yawAngles[0] = 0
            cSet.tiltAngles[0] = angle
            self.assertAlmostEqual(abs(np.radians(angle)), cSet.alphas[0], 8,
                                   'tilt test filled at angle = %f'%angle)
            cSet.tiltAngles[0] = 0

if __name__ == '__main__':
    unittest.main()
