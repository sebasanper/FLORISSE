# -*- coding: utf-8 -*-
import autograd.numpy as np

from turbines.NREL5MW.NREL5MW import NREL5MWTurbine


class Nrel5MWLayout(object):
    """A base layout NREL5MW turbine that will set all turbines to NREL5MW
    when locations are specified. Also sets generic atmospheric conditions"""
    # Atmospheric Conditions
    airDensity = 1.225  # air density
    veer = 0.0          # veer component [deg]
    TI_0 = 0.1          # turbulence intensity [-] ex: 0.1 is 10%
    shear = 0.12        # shear exponent (0.14 -> neutral)

    def __init__(self, usePitch, xLoc, yLoc):
        self._xLoc = xLoc
        self._yLoc = yLoc
        self.nTurbs = len(xLoc)
        self.turbines = [NREL5MWTurbine(usePitch) for i in range(self.nTurbs)]
        self.zLoc = [turb.hubHeight for turb in self.turbines]
        self.windSpeed = 8.0      # wind speed [m/s]
        self.windDirection = 0.0  # wind direction [deg] (compass degrees)

    zLoc = property(lambda self: self._zLoc,
                    lambda self, val: self.setValueAndCheck('_zLoc', val))

    xLoc = property(lambda self: self._xLoc,
                    lambda self, val: self.setValueCheckAndRotate('_xLoc', val))

    yLoc = property(lambda self: self._yLoc,
                    lambda self, val: self.setValueCheckAndRotate('_yLoc', val))

    windDirection = property(lambda self: self._windDirection,
                             lambda self, val: self.setWindDir(val))

    def setWindDir(self, value):
        self._windDirection = value
        self.rotateCoordinates()

    def setValueAndCheck(self, attr, values):
        if len(values) == self.nTurbs:
            setattr(self, attr, np.array(values))
        else:
            raise Exception('length of %s should be %d' % (attr, self.nTurbs))

    def setValueCheckAndRotate(self, attr, values):
        self.setValueAndCheck(attr, values)
        self.rotateCoordinates()


    def rotateCoordinates(self):
        # this function rotates the coordinates so that the flow direction is
        # alligned with the x-axis.
        rotMatrix = self.rotationMatZ(np.radians(self.windDirection))
        coordsRot = np.dot(rotMatrix, np.array([self.xLoc, self.yLoc]))

        # Make the turbine locations start at 0,0. Save the X and Y locations
        self.xLocRot = coordsRot[0, :] - min(coordsRot[0, :])
        self.yLocRot = coordsRot[1, :] - min(coordsRot[1, :])
        self.sortedTurbIds = [i[0] for i in sorted(enumerate(self.xLocRot),
                              key=lambda x:x[1])]

    # Define a rotation matrix around Z
    def rotationMatZ(self, theta):
        return np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]])


class Layout1(Nrel5MWLayout):
    """A windfarm layout with one NREL5MW turbine"""
    def __init__(self, *args):
        super().__init__(*args, [0.0], [0.0])
        # Adjust default atmospheric conditions
        self.windSpeed = 7.0       # wind speed [m/s]


class Layout2by2(Nrel5MWLayout):
    """A windfarm layout with 4 NREL5MW turbines"""
    def __init__(self, *args):
        super().__init__(*args, [0, 800, 0, 800], [0, 0, 600, 600])
        # Adjust default atmospheric conditions
        self.windSpeed = 7.0       # wind speed [m/s]


class Layout3by3(Nrel5MWLayout):
    """A windfarm layout with 9 NREL5MW turbines"""
    def __init__(self, *args):
        xLoc = [300, 300, 300, 1000, 1000, 1000, 1600, 1600, 1600]
        yLoc = [100, 300, 500, 100, 300, 500, 100, 300, 500]
        super().__init__(*args, xLoc, yLoc)


class Layout3(Nrel5MWLayout):
    """A windfarm layout with 3 NREL5MW turbines placed in a line"""
    def __init__(self, *args):
        super().__init__(*args, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
        D = self.turbines[0].rotorDiameter
        self.xLoc = [0, 7*D, 14*D]


class Layout2(Nrel5MWLayout):
    """A windfarm layout with 3 NREL5MW turbines placed in a line"""
    def __init__(self, *args):
        super().__init__(*args, [0.0, 0.0], [0.0, 0.0])
        D = self.turbines[0].rotorDiameter
        self.xLoc = [0, 7*D]


class Layout5(Nrel5MWLayout):
    """A windfarm layout with 5 NREL5MW turbines placed in a line"""
    def __init__(self, *args):
        super().__init__(*args, [0, 630, 1260, 1890, 2520], [0, 0, 0, 0, 0])
