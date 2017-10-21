# -*- coding: utf-8 -*-
# Importing imp causes Spyder to automatically re-import changed modules
import imp
import copy

import florisCoreFunctions.windPlant as windPlant
import optimizers.controlOptimizers as cOpters

import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData

from visualizationTools.viewer import viewer


# Select a velocity, deflection and wake summing model
model = inputClasses.modelData.modelData(2, 1, 2)

# Select a wind farm layout and specify how the turbine control mode
layout = layouts.Layout2(True)
#layout = layouts.Layout2by2(True)
#layout.windDirection = 20
# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.Neutral(layout)

# Run the model and get an output object
outputNeutral = windPlant.windPlant(model, layout, cSet, True)
outputNeutral.printVelocitiesAndPowers()
viewApp = viewer(outputNeutral)
viewApp.showView(0)

# Generate Data to test in model Optimization routine
#powerData = [[[],[],[]],[[],[],[]]]
#layoutPars = {'xLoc': [[0, 500], [0, 800]]}
#cSetPars = {'tiltAngles': ([-20, 0], [0, 0], [20, 0])}
#for keyL, valuesL in layoutPars.items():
#    for i in range(len(valuesL)):
#        for keyC, valuesC in cSetPars.items():
#            for ii in range(len(valuesC)):
#                setattr(layout, keyL, valuesL[i])
#                setattr(cSet, keyC, valuesC[ii])
#                outputNeutral = windPlant.windPlant(model, layout, cSet, True)
#                powerData[i][ii] = outputNeutral.power

#axialOpter = cOpters.axialOptimizer(model, layout, copy.copy(cSet))
#outputOptAxial = axialOpter.optimize()
##viewAppAxialOpt = viewer(outputOptAxial)
##viewAppAxialOpt.showView(0)

#yawOpter = cOpters.yawOptimizer(model, layout, copy.copy(cSet))
#outputOptYaw = yawOpter.optimize()
#viewAppYawOpt = viewer(outputOptYaw)
#viewAppYawOpt.showView(0)

#tiltOpter = cOpters.tiltOptimizer(model, layout, copy.copy(cSet))
#outputOptTilt = tiltOpter.optimize()
#outputOptTilt.cSet.tiltAngles
#outputOptYaw.cSet.yawAngles
##viewAppTiltOpt = viewer(outputOptTilt)
##viewAppTiltOpt.showView(0)

#yawAndAxialOpter = cOpters.yawAndAxialOptimizer(model, layout, copy.copy(cSet))
#outputOptYawAndAxial = yawAndAxialOpter.optimize()
#viewAppYawAndAxialOpt = viewer(outputOptYawAndAxial)
#viewAppYawAndAxialOpt.showView(0)

#from autograd import value_and_grad, grad
#from scipy.optimize import check_grad, approx_fprime
#
#print(check_grad(yawOpter.cost, grad(yawOpter.cost), [0.1, 0.0, 0.0, 0.0, 0.0]))
#print(approx_fprime([0.1, 0.0, 0.0, 0.0, 0.0], yawOpter.cost, 1.49e-08))
#print(value_and_grad(yawOpter.cost)([0.1, 0.0, 0.0, 0.0, 0.0]))
