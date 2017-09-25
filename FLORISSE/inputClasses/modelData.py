#!/usr/bin/python3
# -*- coding: utf-8 -*-
import florisCoreFunctions.wakeCombinationModels as wCMs
import florisCoreFunctions.wakeDeflectionModels as wDMs
import florisCoreFunctions.wakeVelocityModels as wVMs


class modelData:
    """Contains the model parameters that govern wake behaviour"""
    # =========================================================================
    #                               JENSEN and FLORIS Wake Parameters
    # =========================================================================
    wakeDeflection = .17   # standard in literature is 0.17
    wakeExpansion = .05    # wake expansion coefficient
    ad = -4.5              # lateral wake displacement bias parameter (a + bx)
    bd = -0.01             # lateral wake displacement bias parameter (a + bx)
    aT = 0.0               # vertical wake displacement bias parameter (a + bx)
    bT = 0.0               # vertical wake displacement bias parameter
    me = [-0.5, 0.3, 1.0]  # expansion of each region of the wake
    MU = [0.5, 1., 5.5]    # determine velocity of each region in the wake
    aU = 12.0              # wake velocity parameter (a + b*yaw)
    bU = 1.3               # wake velocity parameter (a + b*yaw)

    # =========================================================================
    #                               GAUSS Wake Parameters
    # =========================================================================
    ka = 0.3871            # wake expansion parameter (ka*TI + kb)
    kb = 0.004             # wake expansion parameter (ka*TI + kb)
    alpha = 0.58           # near wake parameter
    beta = 0.077           # near wake parameter

    TIa = 0.73             # magnitude of turbulence added
    TIb = 0.8325           # contribution of turbine operation
    TIc = 0.0325           # contribution of ambient turbulence intensity
    TId = -0.32            # contribution of downstream distance from turbine

    rotorPts = 16
    avgCube = False

    def __init__(self, velocityModel, deflectionModel, combineWakes):
        self.wake = wake

        if velocityModel == 0:
            self.velClass = wVMs.Jensen
        elif velocityModel == 1:
            self.velClass = wVMs.FLORIS
        elif velocityModel == 2:
            self.velClass = wVMs.GAUSS
        else:
            raise NameError('No valid wake velocity model was specified')

        if deflectionModel == 0:
            self.deflClass = wDMs.jimenezDeflection
        elif deflectionModel == 1:
            self.deflClass = wDMs.porteAgelDeflection
        else:
            raise NameError('No valid wake deflection model was specified')

        if combineWakes == 0:
            self.wakeCombine = wCMs.FLS  # freestream linear
        elif combineWakes == 1:
            self.wakeCombine = wCMs.LVLS  # local velocity linear
        elif combineWakes == 2:
            self.wakeCombine = wCMs.SOSFS  # sum of squares freestream
        elif combineWakes == 3:
            self.wakeCombine = wCMs.SOSLVS  # sum of squares local velocity
        else:
            raise NameError('No valid wake combination model was specified')

    # threshold distance of turbines to include in "added turbulence"
    def TIdistance(self, D):
        return 15*D


class wake:
    def __init__(self, model, layout, cSet, output, turbI):
        velocityModel = model.velClass(model, layout, cSet, output, turbI)
        deflectionModel = model.deflClass(model, layout, cSet, output, turbI)
        self.V = velocityModel.V
        self.B = velocityModel.B
        self.displ = deflectionModel.displ
