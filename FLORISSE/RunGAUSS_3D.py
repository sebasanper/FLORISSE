# -*- coding: utf-8 -*-
# Importing imp causes Spyder to automatically re-import changed modules
import imp
import copy

import florisCoreFunctions.windPlant as windPlant
import optimizers.controls.controlOptimizers as optimizers

import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData

from visualizationTools.viewer import viewer


# Select a velocity, deflection and wake summing model
model = inputClasses.modelData.modelData(2, 1, 2)

# Select a wind farm layout and specify how the turbine control mode
layout = layouts.LayoutRow(True)
#layout = layouts.Layout2(True)

# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.Neutral(layout)
#cSet.tiltAngles[0] = -25
#cSet.yawAngles[0] = 25

# Run the model and get an output object
outputNeutral = windPlant.windPlant(model, layout, cSet, True)
outputNeutral.printVelocitiesAndPowers()

#viewApp = viewer(outputNeutral)
#viewApp.showView(0)
#viewApp.showView(3)
#viewApp.showView(4)

# NOTE: large-scale optimization techniques have not been enabled in this
# version, but will be released in future versions
#outputOptim = optimizers.axialOpt(model, layout, copy.copy(cSet))
#outputOptim.printVelocitiesAndPowers()

#outputOptYaw = optimizers.yawOpt(model, layout, copy.copy(cSet))
#outputOptYaw.printVelocitiesAndPowers()
#
#viewAppYawOpt = viewer(outputOptYaw)
#viewAppYawOpt.showView(0)
import autograd.numpy as np
from autograd import value_and_grad

def optPlant(x):
    model.ka = x[0]
    model.kb = x[1]
    output = windPlant.windPlant(model, layout, cSet, False)
    return np.sum(np.array(output.power))

grad_core = value_and_grad(optPlant)
print(grad_core([0.3871, 0.004 ]))
