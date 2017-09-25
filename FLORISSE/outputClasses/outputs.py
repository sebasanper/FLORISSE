# -*- coding: utf-8 -*-
import pickle
import sys
import os.path

from visualizationTools.viewer import viewer


class outputClass:
    """A simple example class"""
    def __init__(self, model, layout, cSet):
        nTurbs = layout.nTurbs
        # Rotated turbines locations
        self.rotLocX = [[] for i in range(nTurbs)]
        self.rotLocY = [[] for i in range(nTurbs)]
        try:
            self.prRootPath = os.path.dirname(sys.modules['__main__'].__file__)
        except:
            self.prRootPath = False

        # individual turbine parameters
        self.Ct = [[] for i in range(nTurbs)]  # Thrust Coefficient
        self.Cp = [[] for i in range(nTurbs)]  # Power Coefficient
        self.aI = [[] for i in range(nTurbs)]  # Axial Induction
        self.TI = [[] for i in range(nTurbs)]  # Turbulence intensity at rotor
        self.windSpeed = [[] for i in range(nTurbs)]  # Windspeed at rotor
        self.power = [[] for i in range(nTurbs)]  # Windspeed at rotor
        self.wakes = [[] for i in range(nTurbs)]  # Windspeed at rotor


class powerOutput(outputClass):
    def __init__(self, model, layout, cSet):
        super().__init__(model, layout, cSet)


class fullOutput(outputClass):
    def __init__(self, model, layout, cSet):
        super().__init__(model, layout, cSet)
        self.model = model
        self.layout = layout
        self.cSet = cSet
        self.viewApp = viewer(self)

    def printVelocitiesAndPowers(self):
        # Velocity and power output
        print('Effective Velocities and turbine powers')
        for i in range(self.layout.nTurbs):
            print('Turbine ', i, ' velocity = ',
                  self.windSpeed[i], ' power = ', self.power[i])
        print(sum(self.power), '\n')
