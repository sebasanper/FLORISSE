#!/usr/bin/python3
# -*- coding: utf-8 -*-


class modelData:
    """A control set for the turbines in the windfarm"""
    # =========================================================================
    #                               FLORIS Wake Parameters
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

    def __init__(self, WakeModel, combineWakes, deflectionModel):
        # 0 = Jensen, 1 = FLORIS, 2 = Gaussian
        self.WakeModel = WakeModel

        # 0 = freestream linear, 1 = local velocity linear
        # 2 = sum of squares freestream, 3 = sum of squares local velocity
        self.combineWakes = combineWakes

        # 0 = Jimenez, 1 = Porte Agel
        self.deflectionModel = deflectionModel
        
    # threshold distance of turbines to include in "added turbulence"
    def TIdistance(self, D):
        return 15*D
