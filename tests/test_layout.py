# -*- coding: utf-8 -*-
import numpy as np
import unittest

import inputClasses.layouts as layouts


class TestLayout(unittest.TestCase):

    """Tests for `FLORISSE\florisCoreFunctions\windPlant.py`."""
    def test_nrel5MWlayout(self):
        layout = layouts.Layout2(True)
        self.assertEqual(layout.airDensity, 1.225)
        self.assertEqual(layout.veer, 0.0)
        self.assertEqual(layout.TI_0, 0.1)
        self.assertEqual(layout.shear, 0.12)
        self.assertEqual(layout.xLoc, [0, 800, 0, 800])
        self.assertEqual(layout.yLoc, [0, 0, 600, 600])
        self.assertEqual(layout.windSpeed, 7.0)
        self.assertEqual(layout.windDirection, 0.0)
        self.assertEqual(layout.nTurbs, len(layout.xLoc))

        for wDir in np.arange(0, 360, 30):
            layout.windDirection = wDir
            xLocRot, yLocRot = self.rotateClockwise(layout.xLoc, layout.yLoc,
                                                    np.radians(wDir))
            self.assertEqual(xLocRot, layout.xLocRot)
            self.assertEqual(yLocRot, layout.yLocRot)

    def rotateClockwise(self, x, y, theta):
        R = np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])
        coordsRot = np.dot(R, np.array([x, y]))
        return (tuple(coordsRot[0, :] - min(coordsRot[0, :])),
                tuple(coordsRot[1, :] - min(coordsRot[1, :])))

if __name__ == '__main__':
    unittest.main()
