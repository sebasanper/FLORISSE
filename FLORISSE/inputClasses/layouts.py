# -*- coding: utf-8 -*-
import numpy as np

from Turbines.NREL5MW.NREL5MW import NREL5MWTurbine


class nrel5MWlayout:
    """A base layout NREL5MW turbine that will set all turbines to NREL5MW
    when locations are specified. Also sets generic atmospheric conditions"""
    # Atmospheric Conditions
    airDensity = 1.225  # air density
    veer = 0.0          # veer component [deg]
    TI_0 = 0.1          # turbulence intensity [-] ex: 0.1 is 10%
    shear = 0.12        # shear exponent (0.14 -> neutral)

    def __init__(self, usePitch):
        self.nTurbs = len(self.xLoc)
        self.turbines = [NREL5MWTurbine(usePitch)
                         for i in range(self.nTurbs)]
        self.zLoc = [turb.hubHeight for turb in self.turbines]

    @property
    def windDirection(self):
        return self.windDir

    # Rotating the wind should change the rotated locations of the turbines.
    @windDirection.setter
    def windDirection(self, windDirection):
        self.windDir = windDirection
        self.rotateCoordinates()

    def rotateCoordinates(self):
        # this function rotates the coordinates so that the flow direction is
        # now alligned with the x-axis. This makes computing wakes and wake
        # overlap much simpler
        rotMatrix = self.rotationMatrixCW(np.radians(self.windDirection))

        coordsRot = np.dot(np.array([self.xLoc, self.yLoc]).T, rotMatrix)

        # Make the turbine locations start at 0,0. Save the X and Y domains
        self.xLocRot = tuple(coordsRot[:, 0].T - min(coordsRot[:, 0]))
        self.yLocRot = tuple(coordsRot[:, 1].T - min(coordsRot[:, 1]))

    # Define a clockwise rotation matrix
    def rotationMatrixCW(self, theta):
        return np.array([[np.cos(theta), np.sin(theta)],
                        [-np.sin(theta),  np.cos(theta)]])


class layout1(nrel5MWlayout):
    """A windfarm layout with one NREL5MW turbine"""
    # set turbine locations - example a single turbine
    xLoc = [0.0]
    yLoc = [0.0]

    def __init__(self, *args):
        super().__init__(*args)
        # Atmospheric Conditions
        self.windSpeed = 7.0            # wind speed [m/s]
        self.windDirection = 0.0       # wind direction [deg] (compass degrees)


class layout2(nrel5MWlayout):
    """A windfarm layout with 4 NREL5MW turbines"""
    # set turbine locations - example 2x2 wind farm
    xLoc = [0, 800, 0, 800]
    yLoc = [0, 0, 600, 600]

    def __init__(self, *args):
        super().__init__(*args)
        # Atmospheric Conditions
        self.windSpeed = 7.0       # wind speed [m/s]
        self.windDirection = 0.0  # wind direction [deg] (compass degrees)


class layout3(nrel5MWlayout):
    """A windfarm layout with 4 NREL5MW turbines"""
    # set turbine locations - example 2x2 wind farm
    xLoc = [300, 300, 300, 1000, 1000, 1000, 1600, 1600, 1600]
    yLoc = [100, 300, 500, 100, 300, 500, 100, 300, 500]

    def __init__(self, *args):
        super().__init__(*args)
        # Atmospheric Conditions
        self.windSpeed = 7.0       # wind speed [m/s]
        self.windDirection = -20.0  # wind direction [deg] (compass degrees)
