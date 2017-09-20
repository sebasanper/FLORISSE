import numpy as np
import scipy 
import matplotlib.pyplot as plt
import pandas as pd
import imp
import os
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import time

import main
import utilities
import wakeModels
import OptModules
import NREL5MW

def determineCpCt(GCPCT,windSpeed):
    
    # general Cp/Ct interpolation function
    
    fCp = interp1d(GCPCT['wind_speed'],GCPCT['CP'])
    fCt = interp1d(GCPCT['wind_speed'],GCPCT['CT'])
    
    Cp = fCp(windSpeed)
    Ct = fCt(windSpeed)
    
    return Cp,Ct

inputData = dict()
# =============================================================================================
#                              Choose Turbine
# =============================================================================================
# turbine info
inputData['TurbineInfo']        = NREL5MW.turbineProperties()

# eventually get rid of these parameters
D     = inputData['TurbineInfo']['rotorDiameter']
HH    = inputData['TurbineInfo']['hubHeight']
beta  = inputData['TurbineInfo']['bladePitch']
gamma = inputData['TurbineInfo']['yawAngle']
phi   = inputData['TurbineInfo']['tilt']
pP    = inputData['TurbineInfo']['pP']
pT    = inputData['TurbineInfo']['pT']
gE    = inputData['TurbineInfo']['generatorEfficiency']
B     = inputData['TurbineInfo']['NumBlades']
TSR   = inputData['TurbineInfo']['TSR']


# =============================================================================================
#                              Turbine Locations
# =============================================================================================
# turbine locations - example 2x2 wind farm
inputData['turbineX'] = [500.0]
inputData['turbineY'] = [250.0]
nTurbs                = len(inputData['turbineX']) 
inputData['turbineZ'] = HH*np.ones(nTurbs)

# =============================================================================================
#                              Atmospheric Conditions
# =============================================================================================
Uinf = 7.0

# atmospheric conditions
inputData['airDensity']          = 1.225           # air density
inputData['windSpeed']           = Uinf            # wind speed [m/s]
inputData['windDirection']       = 270.0           # wind direction [deg] (compass degrees)
inputData['veer']                = 0.0             # veer component [deg]
inputData['turbulenceIntensity'] = 0.1             # turbulence intensity [-] ex: 0.1 is 10% turbulence intensity
inputData['shear']               = 0.12            # shear exponent (0.14 -> neutral)


# =============================================================================================
#                              Turbine Parameters
# =============================================================================================
# individual turbine parameters (update Ct and Cp)
inputData['rotorDiameter']       = [D for i in range(nTurbs)]         # rotor diameter [m]
inputData['hubHeight']           = [HH for i in range(nTurbs)]        # hub height [m]
inputData['pP']                  = [pP for i in range(0,nTurbs)]      # power adjustment to the cosine rule [-]
inputData['pT']                  = [2.07 for i in range(0,nTurbs)]    # power adjustment to the cosine rule [-]
inputData['generatorEfficiency'] = [gE for i in range(nTurbs)]        # generator efficiency [-]
inputData['NumBlades']           = [B for i in range(nTurbs)]         # number of blades [-]
inputData['TSR']                 = [TSR for i in range(nTurbs)]       # tip-speed ratio [-]

#turbine controls
Cp,Ct                           = determineCpCt(inputData['TurbineInfo']['CpCtWindSpeed'],inputData['windSpeed'])
inputData['yawAngles']          = [gamma for i in range(nTurbs)]      # yaw angles [deg]
inputData['tilt']               = [phi for i in range(nTurbs)]        # tilt angle [deg]
inputData['Cp']                 = [Cp for i in range(nTurbs)]         # power coefficient [-]
inputData['bladePitch']         = [beta for i in range(nTurbs)]       # collective blade pitch angle 
inputData['Ct']                 = [Ct for i in range(nTurbs)]         # thrust coefficient [-] 

# =============================================================================================
#                               Choose Wake Model
# =============================================================================================
# wake models
inputData['WakeModel']          = 2                                   # 0 = Jensen, 1 = FLORIS, 2 = Gaussian

# =============================================================================================
#                               Choose wake combination
# =============================================================================================
# velocity deficity
inputData['combineWakes']       = 2                                   # 0 = freestream linear, 1 = local velocity linear, 2 = sum of squares freestream, 3 = sum of squares local velocity 


# =============================================================================================
#                               Wake Parameters (FLORIS, GAUSS)
# =============================================================================================
# wake parameters (Jensen and FLORIS)
inputData['wakeDeflection'] = [0.17 for i in range(nTurbs)]    # standard in literature is 0.17
inputData['wakeExpansion']  = [0.05 for i in range(nTurbs)]    # wake expansion coefficient
inputData['ad']             = -4.5                             # lateral wake displacement bias parameter (a + bx)
inputData['bd']             = -0.01                            # lateral wake displacement bias parameter (a + bx)
inputData['aT']             = 0.0                              # vertical wake displacement bias parameter (a + bx)
inputData['bT']             = 0.0                              # vertical wake displacement bias parameter
inputData['me']             = [-0.5,0.3,1.0]                   # expansion of each region of the wake
inputData['MU']             = [0.5,1.,5.5]                     # determine velocity of each region in the wake
inputData['aU']             = 12.0                             # wake velocity parameter (a + b*yaw)                                                    
inputData['bU']             = 1.3                              # wake velocity parameter (a + b*yaw) 

# wake parameters GAUSS
inputData['ka'] = 0.3871                    # wake expansion parameter (ka*TI + kb)                     
inputData['kb'] = 0.004                     # wake expansion parameter (ka*TI + kb)
inputData['alpha'] = 0.58                   # near wake parameter
inputData['beta'] = 0.077                   # near wake parameter

inputData['TIdistance'] = 15*D              # threshold distance of turbines to include in "added turbulence"
inputData['TIa'] = 0.73                     # magnitude of turbulence added
inputData['TIb'] = 0.8325                   # contribution of turbine operation
inputData['TIc'] = 0.0325                   # contribution of ambient turbulence intensity
inputData['TId'] = -0.32                    # contribution of downstream distance from turbine

# =============================================================================================
#                                      Optimization
# =============================================================================================

# optimization options
inputData['axial_opt']        = False       # True turns on thrust control
inputData['yaw_opt']          = False       # True turns on wake steering

# NOTE: large-scale optimization techniques have not been enabled in this version, but will be released in future versions

# constraints on controls
inputData['minYaw']           = [0.0]        # minimum yaw angles considered in wake steering optimization
inputData['maxYaw']           = [25.0]       # maximum yaw angles considered in wake steering optimization

# =============================================================================================
#                                       Visualization
# =============================================================================================
# domain for plotting
inputData['rotorPts']           = 16          # number of points evaluated on the turbine rotor 
inputData['avgCube']            = False       # True: take cube root of the average of cubed wind speed, False: take the average of the wind speed

inputData['xLen']               = [np.min(inputData['turbineX'])-2*D,np.max(inputData['turbineX'])+15*D] # x domain (min,max)
inputData['yLen']               = [np.min(inputData['turbineY'])-2*D,np.max(inputData['turbineY'])+2*D]  # y domain (min,max)
inputData['zLen'] 				= [0.0,2*HH]                                                             # z domain (min,max)

inputData['outputFlowField']     = False        # output full flow field
inputData['visualizeHorizontal'] = False        # visualize horizontal flow field at hub height       
inputData['nSamplesX']           = 200          # resolution of the plot in the x direction
inputData['nSamplesY']           = 200          # resolution of the plot in the y direction
inputData['nSamplesZ']           = 50           # resoltuion of the plot in the z direction

# plot cut through slices (parallel with the rotor)
inputData['cutThrough']         = False
inputData['cutTurbID']          = 0
inputData['downLocs']           = [2*D]

# =============================================================================================
#                                       Lidar
# =============================================================================================

# Lidar module based on the first turbine - University of Stuttgard Lidar model
# only one lidar can be used currently for one turbine
inputData['Lidar']              = False                         # Turn on lidar module and output flow field from 5 downstream locations 
inputData['turbineLidar']       = 0                             # identify the turbine that the Lidar is on
inputData['xLidar']             = pd.read_csv('LiDAR/xLoc.csv') # x locations of the lidar measurements
inputData['yLidar']             = pd.read_csv('LiDAR/yLoc.csv') # y locations of the lidar measurements
inputData['zLidar']             = pd.read_csv('LiDAR/zLoc.csv') # z locations of the lidar measurements
inputData['LidarParams']        = 'LiDAR/Liss2Grid7x7widescreen_1d000_0d450_2d800_10000_1d600_NRELGE_Continuous_Parameter_JENExport.mat'
inputData['LidarTrajectory']    = 'LiDAR/Liss2Grid7x7widescreen_1d000_0d450_2d800_10000_1d600_NRELGE_Continuous_JENExport.mat'

# =============================================================================================
#                                      Output points
# =============================================================================================
# select points to output velocity (if you want the full flow field, you should set outputFlowField to True rather than use this)
inputData['xPts'] = np.concatenate(inputData['xLidar'].values) + inputData['turbineX'][inputData['turbineLidar']]
inputData['yPts'] = np.concatenate(inputData['yLidar'].values) + inputData['turbineY'][inputData['turbineLidar']]
inputData['zPts'] = np.concatenate(inputData['zLidar'].values) + inputData['turbineZ'][inputData['turbineLidar']]
inputData['points'] = False    # must set this to true if you want points


# =============================================================================================
#                                      example (basic)
# =============================================================================================
# get velocity and power out
inputData['visualizeHorizontal'] = False
inputData['yaw_opt']             = False
inputData['axial_opt']           = False
inputData['Lidar']               = False
inputData['points']              = False
outputData = main.windPlant(inputData)

print('Effective Velocities')
for i in range(nTurbs):
    print('Turbine ', i, ' velocity = ', outputData['Ueff'][i])

print('Power')
for i in range(nTurbs):
    print('Turbine ', i, ' power = ', outputData['powerOut'][i])

# =============================================================================================
#                          example (visualization Porte-Agel wake)
# =============================================================================================
# Cut a slice, unyawed and yawed conditions
imp.reload(main)

inputData['visualizeHorizontal'] = True
inputData['yaw_opt']             = False
inputData['axial_opt']           = False
inputData['Lidar']               = False
inputData['points']              = False
inputData['WakeModel'] = 2
inputData['yawAngles'] = np.zeros(nTurbs)
#outputData = main.windPlant(inputData)

inputData['visualizeHorizontal'] = True
inputData['yawAngles'][0] = 25.0
outputData = main.windPlant(inputData)

# =============================================================================================
#                          example (visualization Zoned wake)
# =============================================================================================
# Cut a slice, unyawed and yawed conditions
inputData['visualizeHorizontal'] = True
inputData['yaw_opt']             = False
inputData['axial_opt']           = False
inputData['Lidar']               = False
inputData['points']              = False
inputData['WakeModel'] = 1
inputData['yawAngles'] = np.zeros(nTurbs)
outputData = main.windPlant(inputData)

inputData['visualizeHorizontal'] = True
inputData['yawAngles'][0] = 25.0
outputData = main.windPlant(inputData)
