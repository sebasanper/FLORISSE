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
        self.rootPath = os.path.dirname(sys.modules['__main__'].__file__)

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

    def writeToDisk(self, filename):
        vDataPath = self.rootPath + '/visualizationData/'
        fullFileName = vDataPath+filename
        os.rename(self.v.vtiFile, fullFileName+'.vti')
        self.v.vtiFile = fullFileName+'.vti'
        pickle.dump(self, open(fullFileName+'.p', "wb" ))
