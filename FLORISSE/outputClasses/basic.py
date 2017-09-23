# -*- coding: utf-8 -*-

import numpy as np


class output:
    """A simple example class"""
    def __init__(self, layout):
        nTurbs = layout.nTurbs

        # individual turbine parameters (update Ct and Cp)
        self.Ct = [0 for i in range(nTurbs)]  # Thrust Coefficient
        self.Cp = [0 for i in range(nTurbs)]  # Power Coefficient
        self.aI = [0 for i in range(nTurbs)]  # Axial Induction
        self.TI = [0 for i in range(nTurbs)]  # Turbulence intensity at rotor
        self.windSpeed = [0 for i in range(nTurbs)]  # Windspeed at rotor
        self.power = [0 for i in range(nTurbs)]  # Windspeed at rotor

    def computePower(self, layout, cSet):
        yaw = cSet.yawAngles
        tilt = cSet.tiltAngles

        for i in range(layout.nTurbs):
            turb = layout.turbines[i]
            Cptmp = (self.Cp[i]*(np.cos(yaw[i]*np.pi/180.)**turb.pP) *
                     (np.cos((tilt[i])*np.pi/180.)**turb.pT))
            self.power[i] = (.5*layout.airDensity*(np.pi*(turb.rotorDiameter/2)**2) *
                             Cptmp*turb.generatorEfficiency*self.windSpeed[i]**3)
