# -*- coding: utf-8 -*-
import numpy as np
import unittest

import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData


class TestControlSet(unittest.TestCase):

    """Tests for `FLORISSE\florisCoreFunctions\windPlant.py`."""
    def test_neutral_controlSet(self):
        layout = layouts.Layout1(True)
        cSet = inputClasses.controlSettings.Neutral(layout)

        with self.assertRaises(Exception):
            cSet.yawAngles = [0, 0, 0, 0, 0]
        with self.assertRaises(IndexError):
            cSet.yawAngles[layout.nTurbs+1]

        for i in range(layout.nTurbs):
            self.assertEqual(cSet.yawAngles[i], 0)
            self.assertEqual(cSet.tiltAngles[i], 0)

        for angle in np.arange(-30, 30, 15):
            # TODO: don't use abs() but sign derived from wakeDir
            cSet.yawAngles[0] = angle
            self.assertAlmostEqual(abs(np.radians(angle)), cSet.phis[0], 8,
                                   'yaw test filled at angle = %f' % angle)
            cSet.yawAngles[0] = 0
            cSet.tiltAngles[0] = angle
            self.assertAlmostEqual(abs(np.radians(angle)), cSet.phis[0], 8,
                                   'tilt test filled at angle = %f' % angle)
            cSet.tiltAngles[0] = 0

    def test_deflection_direction(self):
        layout = layouts.Layout1(True)
        cSet = inputClasses.controlSettings.Neutral(layout)
        # Generate and test a wake
        for yaw in [-20, 0, 20]:
            for tilt in [-20, 0, 20]:
                cSet.yawAngles[0] = yaw
                cSet.tiltAngles[0] = tilt
                self.assertEqual(np.sign(yaw), -np.sign(cSet.wakeDir[0][1]))
                self.assertEqual(np.sign(tilt), -np.sign(cSet.wakeDir[0][2]))

if __name__ == '__main__':
    unittest.main()
