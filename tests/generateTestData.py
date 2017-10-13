# -*- coding: utf-8 -*-
import imp
import pickle
import autograd.numpy as np

import florisCoreFunctions.windPlant as windPlant

import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData

from turbines.NREL5MW.NREL5MW import NREL5MWTurbine


# Select a velocity, deflection and wake summing model
model = inputClasses.modelData.modelData(2, 1, 2)

# Select a wind farm layout and specify how the turbine control mode
layout = layouts.Layout2(True)
turbine = NREL5MWTurbine(True)

# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.Neutral(layout)

powerDataWindRose2by2 = {}
for deg in range(0, 360, 20):
    layout.windDirection = deg
    outputObject = windPlant.windPlant(model, layout, cSet, False)
    powerDataWindRose2by2[deg] = sum(outputObject.power)[0]

pickle.dump(powerDataWindRose2by2, open('powerDataWindRose2by2.p', 'wb'))

CpCtWithPitch = {}
for pitch in np.arange(turbine.betaLims[0], turbine.betaLims[1], 1):
    CpCtWithPitch[pitch] = {}
    for windSpeed in np.arange(0, 20, 1):
        CpCtWithPitch[pitch][windSpeed] = turbine.Ct(windSpeed, pitch)
pickle.dump(CpCtWithPitch, open('CpCtWithPitch.p', 'wb'))

turbine = NREL5MWTurbine(False)
CpCtWithoutPitch = {}
for windSpeed in np.arange(0, 20, .1):
    CpCtWithoutPitch[windSpeed] = turbine.Ct(windSpeed)
pickle.dump(CpCtWithoutPitch, open('CpCtWithoutPitch.p', 'wb'))

