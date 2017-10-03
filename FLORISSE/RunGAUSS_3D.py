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
layout = layouts.Layout2(True)

# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.Neutral(layout)
#cSet.tiltAngles[0] = -25
#cSet.yawAngles[0] = 25

# Run the model and get an output object
outputNeutral = windPlant.windPlant(model, layout, cSet, True)
outputNeutral.printVelocitiesAndPowers()

viewApp = viewer(outputNeutral)
viewApp.showView(0)
#viewApp.showView(3)
#viewApp.showView(4)

# NOTE: large-scale optimization techniques have not been enabled in this
# version, but will be released in future versions
#outputOptim = optimizers.axialOpt(model, layout, copy.copy(cSet))
#outputOptim.printVelocitiesAndPowers()

outputOptYaw = optimizers.yawOpt(model, layout, copy.copy(cSet))
outputOptYaw.printVelocitiesAndPowers()

viewAppYawOpt = viewer(outputOptYaw)
viewAppYawOpt.showView(0)

#import numpy as np
#a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
#print(np.einsum('ij->i', a))
#print(np.einsum('ji->i', a))
#print(np.einsum('ij->j', a))
#print(np.einsum('ji->j', a))