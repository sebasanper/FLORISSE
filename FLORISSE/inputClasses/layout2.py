#!/usr/bin/python3
# -*- coding: utf-8 -*-

from Turbines.NREL5MW.NREL5MW import NREL5MWTurbine


class layout:
    """A windfarm layout with 4 NREL5MW turbines"""
    # set turbine locations - example 2x2 wind farm
    turbineX = [0, 500, 0, 500]
    turbineY = [0, 0, 500, 500]

    # Atmospheric Conditions
    windSpeed = 7.0            # wind speed [m/s]
    airDensity = 1.225         # air density
    windDirection = 0.0       # wind direction [deg] (compass degrees)
    veer = 0.0                 # veer component [deg]
    turbulenceIntensity = 0.1  # turbulence intensity [-] ex: 0.1 is 10%
    shear = 0.12               # shear exponent (0.14 -> neutral)

    def __init__(self, usePitch, useTSR):
        self.nTurbs = len(self.turbineX)
        self.turbines = [NREL5MWTurbine(usePitch, useTSR)
                         for i in range(self.nTurbs)]
        self.turbineZ = [turb.hubHeight
                         for turb in self.turbines]
