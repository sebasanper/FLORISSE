#!/usr/bin/python3
# -*- coding: utf-8 -*-

from Turbines.NREL5MW.NREL5MW import NREL5MWTurbine


class nrel5MWlayout:
    """A base layout NREL5MW turbine that will set all turbines to NREL5MW
    when locations are specified. Also sets generic atmospheric conditions"""
    # Atmospheric Conditions
    airDensity = 1.225  # air density
    veer = 0.0          # veer component [deg]
    TI_0 = 0.1          # turbulence intensity [-] ex: 0.1 is 10%
    shear = 0.12        # shear exponent (0.14 -> neutral)

    def __init__(self, usePitch, useTSR):
        self.nTurbs = len(self.locX)
        self.turbines = [NREL5MWTurbine(usePitch, useTSR)
                         for i in range(self.nTurbs)]
        self.locZ = [turb.hubHeight for turb in self.turbines]


class layout1(nrel5MWlayout):
    """A windfarm layout with one NREL5MW turbine"""
    # set turbine locations - example a single turbine
    locX = [0.0]
    locY = [0.0]

    # Atmospheric Conditions
    windSpeed = 7.0            # wind speed [m/s]
    windDirection = 0.0       # wind direction [deg] (compass degrees)

    def __init__(self, usePitch, useTSR):
        super().__init__(usePitch, useTSR)


class layout2(nrel5MWlayout):
    """A windfarm layout with 4 NREL5MW turbines"""
    # set turbine locations - example 2x2 wind farm
    locX = [0, 800, 0, 800]
    locY = [0, 0, 600, 600]

    # Atmospheric Conditions
    windSpeed = 7.0            # wind speed [m/s]
    windDirection = 0       # wind direction [deg] (compass degrees)

    def __init__(self, usePitch, useTSR):
        super().__init__(usePitch, useTSR)


class layout3(nrel5MWlayout):
    """A windfarm layout with 4 NREL5MW turbines"""
    # set turbine locations - example 2x2 wind farm
    locX = [0, 500, 1000, 0, 500, 1000, 0, 500, 1000]
    locY = [0, 0, 0, 500, 500, 500, 1000, 1000, 1000]

    # Atmospheric Conditions
    windSpeed = 7.0            # wind speed [m/s]
    windDirection = 0       # wind direction [deg] (compass degrees)

    def __init__(self, usePitch, useTSR):
        super().__init__(usePitch, useTSR)
