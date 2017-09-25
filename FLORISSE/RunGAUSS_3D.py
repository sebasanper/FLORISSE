# Importing imp causes Spyder to automatically re-import changed modules
import imp
import copy

import florisCoreFunctions.windPlant as windPlant
import florisCoreFunctions.OptModules as optimizers

from inputClasses.layouts import layout2 as layoutClass
import inputClasses.controlSettings
import inputClasses.modelData


# Select a velocity, deflection and wake summing model
model = inputClasses.modelData.modelData(2, 1, 2)

# Select a wind farm layout and specify how the turbine control mode
layout = layoutClass(True, False)

# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.neutral(layout)

# Run the model and get an output object
outputNeutral = windPlant.windPlant(model, layout, cSet, True)
outputNeutral.printVelocitiesAndPowers()
#outputNeutral.viewApp.showView(0)

# NOTE: large-scale optimization techniques have not been enabled in this
# version, but will be released in future versions
#outputOptim = optimizers.axialOpt(model, layout, copy.copy(cSet))
#outputOptim.printVelocitiesAndPowers()
#outputOptim.viewApp.showView(0)
#
#outputOptim = optimizers.yawOpt(model, layout, copy.copy(cSet))
#outputOptim.printVelocitiesAndPowers()
#outputOptim.viewApp.showView(0)
