# -*- coding: utf-8 -*-
import imp
import copy

import autograd.numpy as np

import florisCoreFunctions.windPlant as windPlant
import optimizers.modelOptimizers as mOpters
import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData
#import optimizationData.tiltHeigthPowerData as tiltHeigthPowerData

# Select a velocity, deflection and wake summing model
model = inputClasses.modelData.modelData(2, 1, 2)

# Select a wind farm layout and specify how the turbine control mode
layout = layouts.Layout2(True)

# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.Neutral(layout)
outputNeutral = windPlant.windPlant(model, layout, cSet, True)
neutralPower = sum(outputNeutral.power)

#powerTargets = tiltHeigthPowerData.powerTargets
#layoutPars = tiltHeigthPowerData.layoutPars
#cSetPars = tiltHeigthPowerData.cSetPars
powerTargets = (([1525722, 817804],
                 [1735380, 544703],
                 [1525722, 834347]),
                ([1525722, 1191924],
                 [1735380, 920208],
                 [1525722, 1201671]))
layoutPars = {'xLoc': [[0, 500], [0, 800]]}
cSetPars = {'tiltAngles': ([-20, 0], [0, 0], [20, 0])}

alphaBetaOpter = mOpters.alphaBetaOptimizer(copy.copy(model), layout, cSet,
                                            powerTargets, layoutPars, cSetPars)
outputOptAlphaBeta = alphaBetaOpter.optimize()

print(outputOptAlphaBeta.power)
print([outputOptAlphaBeta.model.alpha, outputOptAlphaBeta.model.beta])
print([model.alpha, model.beta], '\n')

kaKbOpter = mOpters.kaKbOptimizer(copy.copy(model), layout, cSet,
                                  powerTargets, layoutPars, cSetPars)
outputOptKaKb = kaKbOpter.optimize()

print(outputOptKaKb.power)
print([outputOptKaKb.model.ka.value, outputOptKaKb.model.kb.value,
       outputOptKaKb.model.ad.value, outputOptKaKb.model.bd.value,
       outputOptKaKb.model.aT.value, outputOptKaKb.model.bT.value])
print([model.ka, model.kb, model.ad, model.bd, model.aT, model.bT], '\n')

TIOpter = mOpters.TIOptimizer(copy.copy(model), layout, cSet,
                              powerTargets, layoutPars, cSetPars)
outputOptTI = TIOpter.optimize()

print(outputOptTI.power)
print([outputOptTI.model.TIa.value, outputOptTI.model.TIb.value,
       outputOptTI.model.TIc.value, outputOptTI.model.TId.value])
print([model.TIa, model.TIb, model.TIc, model.TId], '\n')
