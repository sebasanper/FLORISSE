# -*- coding: utf-8 -*-
import sys
import os.path


class outputClass:
    """A simple example class"""
    def __init__(self, model, layout, cSet):

        # Try to extract the project path so flowfield data can be stored
        try:
            if hasattr(sys.modules['__main__'], '__file__'):
                self.prPath = os.path.dirname(sys.modules['__main__'].__file__)
            else:
                curPath = sys.modules['outputClasses.outputs'].__file__
                pathAsList = os.path.dirname(curPath).split(os.sep)[:-1]
                self.prPath = os.sep.join(pathAsList)
        except:
            self.prPath = False

        # individual turbine parameters
        self.Ct = []  # Thrust Coefficient
        self.Cp = []  # Power Coefficient
        self.aI = []  # Axial Induction
        self.TI = []  # Turbulence intensity at rotor
        self.windSpeed = []  # Windspeed at rotor
        self.power = []  # Windspeed at rotor
        self.wakes = []  # Windspeed at rotor


class powerOutput(outputClass):
    def __init__(self, model, layout, cSet):
        super().__init__(model, layout, cSet)


class fullOutput(outputClass):
    def __init__(self, model, layout, cSet):
        super().__init__(model, layout, cSet)
        self.model = model
        self.layout = layout
        self.cSet = cSet

    def printVelocitiesAndPowers(self):
        # Velocity and power output
        print('Effective Velocities and turbine powers')
        for i in range(self.layout.nTurbs):
            print('Turbine ', i, ' velocity = ',
                  self.windSpeed[i], ' power = ', self.power[i])
        print(sum(self.power), '\n')
