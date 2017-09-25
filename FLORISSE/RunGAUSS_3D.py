# Importing imp causes Spyder to automatically re-import changed modules
import imp

import florisCoreFunctions.windPlant as windPlant
from inputClasses.layouts import layout2 as layoutClass
import inputClasses.controlSettings
import inputClasses.modelData


# =============================================================================
#                      Basic model run for power prediction
# =============================================================================

# Select a velocity, deflection and wake combining model
model = inputClasses.modelData.modelData(2, 1, 2)

# Select a wind farm layout and specify how the turbines in that layout behave
layout = layoutClass(True, False)

# Generate control settings for the turbines in the layout
# all turbines set aligned with wind
cSet = inputClasses.controlSettings.neutral(layout)

# Run the model and get an output object
output = windPlant.windPlant(model, layout, cSet, True)

# Power output
for i in range(layout.nTurbs):
    print(output.windSpeed[i])
    print(output.power[i])
print(sum(output.power))

output.viewApp.showView(3)

#import pickle
#outputOldOne = pickle.load( open('C:/Users/roald/Dropbox/Afstuderen/'+
#          'Python/FLORISSE/FLORISSE/visualizationData/simpleRun.p', "rb" ))
#outputOldOne.v.showView(4)

# =============================================================================================
#                                       Visualization
# =============================================================================================
#vis = dict()
#D = layout.TurbineInfo['rotorDiameter']

#outputView = viewer(output)
#outputView.showView(4)
#file_Name = "testfile"
#
#pickle.dump( favorite_color, open( "save.p", "wb" ) )
#favorite_color = pickle.load( open( "save.p", "rb" ) )

#visualizeOutput(output)

# plot cut through slices (parallel with the rotor)
#vis['cutThrough'] = False
#vis['cutTurbID'] = 0
#vis['downLocs'] = [2*D]
#
## output full flow field
#vis['outputFlowField'] = False
## visualize horizontal flow field at hub height
#vis['visualizeHorizontal'] = False
#
#vis['xLen'] = [min(layout.turbineX)-2*D, max(layout.turbineX)+15*D]
#vis['yLen'] = [min(layout.turbineY)-2*D, max(layout.turbineY)+2*D]
#vis['zLen'] = [0.0, 2*max(layout.turbineZ)]
#
#vis['nSamplesX'] = 200          # resolution of the plot in the x direction
#vis['nSamplesY'] = 200          # resolution of the plot in the y direction
#vis['nSamplesZ'] = 50           # resoltuion of the plot in the z direction
#
#
## Lidar module based on the first turbine - University of Stuttgard Lidar model
## only one lidar can be used currently for one turbine
#
## lidar module plots flow field at 5 downstream locations
#vis['Lidar'] = False
#vis['lidarTubrID'] = 0           # identify the turbine that the Lidar is on
#vis['xLidar'] = read_csv('LiDAR/xLoc.csv') # x locations of the lidar measurements
#vis['yLidar'] = read_csv('LiDAR/yLoc.csv') # y locations of the lidar measurements
#vis['zLidar'] = read_csv('LiDAR/zLoc.csv') # z locations of the lidar measurements
#vis['LidarParams'] = 'LiDAR/Liss2Grid7x7widescreen_1d000_0d450_2d800_10000_1d600_NRELGE_Continuous_Parameter_JENExport.mat'
#vis['LidarTrajectory'] = 'LiDAR/Liss2Grid7x7widescreen_1d000_0d450_2d800_10000_1d600_NRELGE_Continuous_JENExport.mat'
#
### select points to output velocity (if you want the full flow field, you should set outputFlowField to True rather than use this)
##inputData['xPts'] = np.concatenate(inputData['xLidar'].values) + inputData['turbineX'][inputData['turbineLidar']]
##inputData['yPts'] = np.concatenate(inputData['yLidar'].values) + inputData['turbineY'][inputData['turbineLidar']]
##inputData['zPts'] = np.concatenate(inputData['zLidar'].values) + inputData['turbineZ'][inputData['turbineLidar']]
##inputData['points'] = False    # must set this to true if you want points

## =============================================================================================
##                                      Optimization
## =============================================================================================
## NOTE: large-scale optimization techniques have not been enabled in this version, but will be released in future versions
#optim = dict()
#
## optimization options
#optim['axial_opt'] = False  # True turns on thrust control
#optim['yaw_opt'] = False    # True turns on wake steering
#

## =============================================================================================
##                                      example (basic)
## =============================================================================================
## get velocity and power out
#
#outputData = main.windPlant(model, layout, cSet)
#
## plot the flow field if specified
#if inputData['visualizeHorizontal']:
#    MPLVisualizations.visualizeHorizontal(xLen,yLen,zLen,Ufield,inputData)
#elif inputData['cutThrough']:
#    MPLVisualizations.visualizeCut(xLen,yLen,zLen,Ufield,inputData)
#elif inputData['Lidar']:
#    Upts = utilities.outputUpts(inputData,X,Y,Z,Ufield)
#    vlos = utilities.VLOS(x_W,y_W,z_W,inputData,Upts)
#    #MPLVisualizations.visualizeLidar(xTurb[idxLidar],yTurb[idxLidar],X,Y,Z,Ufield,inputData,vlos)
#
#if inputData['outputFlowField']:
#    return Ueff,Ufield,xLen,yLen
#elif inputData['Lidar']:
#    return Ueff,Ufield,X,Y,Z,Upts,vlos
#elif inputData['points']:
#    Upts = utilities.outputUpts(inputData,X,Y,Z,Ufield)
#    return Ueff, Upts
#else:
#    return Ueff
#    
#print('Effective Velocities')
#for i in range(nTurbs):
#    print('Turbine ', i, ' velocity = ', outputData['Ueff'][i])
#
#print('Power')
#for i in range(nTurbs):
#    print('Turbine ', i, ' power = ', outputData['powerOut'][i])
#
## =============================================================================================
##                          Horizontal visualization Porte-Agel wake
## =============================================================================================
## Cut a slice, unyawed and yawed conditions
#
#inputData['visualizeHorizontal'] = True
#inputData['yaw_opt']             = False
#inputData['axial_opt']           = False
#inputData['Lidar']               = False
#inputData['points']              = False
#inputData['WakeModel'] = 2
#inputData['yawAngles'] = np.zeros(nTurbs)
##outputData = main.windPlant(inputData)
#
#inputData['visualizeHorizontal'] = True
#inputData['yawAngles'][0] = 25.0
#outputData = main.windPlant(inputData)
#
## =============================================================================================
##                          Horizontal visualization Zoned wake
## =============================================================================================
## Cut a slice, unyawed and yawed conditions
##inputData['visualizeHorizontal'] = True
##inputData['yaw_opt']             = False
##inputData['axial_opt']           = False
##inputData['Lidar']               = False
##inputData['points']              = False
##inputData['WakeModel'] = 1
##inputData['yawAngles'] = np.zeros(nTurbs)
##outputData = main.windPlant(inputData)
##
##inputData['visualizeHorizontal'] = True
##inputData['yawAngles'][0] = 25.0
##outputData = main.windPlant(inputData)
#
## =============================================================================================
##                          Lidar visualization Porte-Agel wake
## =============================================================================================
## Cut a slice, unyawed and yawed conditions
#
#inputData['visualizeHorizontal'] = false
#inputData['yaw_opt']             = False
#inputData['axial_opt']           = False
#inputData['Lidar']               = False
#inputData['points']              = False
#inputData['WakeModel'] = 2
#inputData['yawAngles'] = np.zeros(nTurbs)
##outputData = main.windPlant(inputData)
#
#inputData['visualizeHorizontal'] = True
#inputData['yawAngles'][0] = 25.0
#outputData = main.windPlant(inputData)