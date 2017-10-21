# -*- coding: utf-8 -*-
import numpy as np
import unittest

import inputClasses.layouts as layouts


class TestLayout(unittest.TestCase):

    """Tests for `FLORISSE\florisCoreFunctions\windPlant.py`."""
    def test_nrel5MWlayout(self):
        layout = layouts.Layout2by2(True)
        self.assertEqual(layout.airDensity, 1.225)
        self.assertEqual(layout.veer, 0.0)
        self.assertEqual(layout.TI_0, 0.1)
        self.assertEqual(layout.shear, 0.12)
        np.testing.assert_array_equal(layout.xLoc, np.array([0, 800, 0, 800]))
        np.testing.assert_array_equal(layout.yLoc, np.array([0, 0, 600, 600]))
        self.assertEqual(layout.windSpeed, 7.0)
        self.assertEqual(layout.windDirection, 0.0)
        self.assertEqual(layout.nTurbs, len(layout.xLoc))

        # TODO: change this test to test some usefell edge cases instead of a
        # vector with a bunch of values
        for wDir in np.arange(0, 360, 30):
            layout.windDirection = wDir
            xLocRot, yLocRot = self.rotateClockwise(layout.xLoc, layout.yLoc,
                                                    np.radians(wDir))
            np.testing.assert_array_equal(xLocRot, layout.xLocRot)
            np.testing.assert_array_equal(yLocRot, layout.yLocRot)

    def rotateClockwise(self, x, y, theta):
        R = np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])
        coordsRot = np.dot(R, np.array([x, y]))
        return (coordsRot[0, :] - min(coordsRot[0, :]),
                coordsRot[1, :] - min(coordsRot[1, :]))

if __name__ == '__main__':
    unittest.main()
