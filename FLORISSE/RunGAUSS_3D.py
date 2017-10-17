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
#layout = layouts.LayoutRow(True)
layout = layouts.Layout2(True)

# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.Neutral(layout)

# Run the model and get an output object
#outputNeutral = windPlant.windPlant(model, layout, cSet, True)
#outputNeutral.printVelocitiesAndPowers()

#viewApp = viewer(outputNeutral)
#viewApp.showView(0)
#viewApp.showView(3)
#viewApp.showView(4)

#outputOptYaw = optimizers.yawOpt(model, layout, copy.copy(cSet))
#outputOptYaw.printVelocitiesAndPowers()

#viewAppYawOpt = viewer(outputOptYaw)
#viewAppYawOpt.showView(0)
from autograd import value_and_grad


def KaKbCost(x):
    model.alpha = x[0]
    model.beta = x[1]
    output = windPlant.windPlant(model, layout, cSet, False)
    return sum(output.power)

print(value_and_grad(KaKbCost)([0.3871, 0.004]))


def AlphaBetaCost(x):
    model.alpha = x[0]
    model.beta = x[1]
    output = windPlant.windPlant(model, layout, cSet, False)
    return sum(output.power)

print(value_and_grad(AlphaBetaCost)([0.58, 0.077]))


def TICost(x):
    model.TIa = x[0]
    model.TIb = x[1]
    model.TIc = x[2]
    model.TId = x[3]
    output = windPlant.windPlant(model, layout, cSet, False)
    return sum(output.power)

print(value_and_grad(TICost)([0.73, 0.8325, 0.0325, -0.32]))


def yawCost(x):
    cSet.yawAngles = x
    output = windPlant.windPlant(model, layout, cSet, False)
    return sum(output.power)

print(value_and_grad(yawCost)([0.1, 0.0, 0.0, 0.0]))
