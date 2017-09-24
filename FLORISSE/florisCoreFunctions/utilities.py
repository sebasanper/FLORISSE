# ==============================================================================
# Utilities for wake model
# ==============================================================================

import numpy as np
import copy
from scipy.interpolate import griddata
import scipy.io
import florisCoreFunctions.wakeModels as wakeModels

# ==================================================================FUNCTIONS============================================================================
# 1 getWindCoords
# 2 getPointsVLOS - retrieves the line of sight wind speed (vlos) from a wind field
# 3 VLOS
# 4 rotatedCoordinates - rotate the turbines to 270 degrees and evaluate the wake model in this way.  Makes coding the wake models simpler
# 5 avgVelocity - compute the average velocity across a rotor
# 6 combineWakes - combine the wakes based linear superposition or sum of squares
# 7 computeAddedTI - compute the turbulence intensity added to the wake based on turbine operation and ambient turbulence intensity
# 8 computeOverlap - compute overlap of an upstream turbine wake (based on a specified threshold)
# 9 computePower - compute the power at each turbine based on the effective velocity determined using the wake model
# 10 initializeFlowField - initialize the flow field used in the 3D model based on shear using the power log law
# 11 determineCpCt - interpolation function based on Cp, Ct, and wind speed
# 12 outputUpts - compute the output velocity at the points specified in inputData
# ======================================================================================================================================================


def getWindCoords(inputData):

    Parameter = scipy.io.loadmat(inputData['LidarParams'])
    Trajectory = scipy.io.loadmat(inputData['LidarTrajectory']) 

    a = Parameter['Parameter'][0][0][0]['a'][0][0][0]
    Yaw = Parameter['Parameter'][0][0][0]['YawAngle'][0][0][0]
    Pitch = Parameter['Parameter'][0][0][0]['PitchAngle'][0][0][0]
    Roll = Parameter['Parameter'][0][0][0]['RollAngle'][0][0][0]
    PositionLinW = Parameter['Parameter'][0][0][0]['PositionLinI'][0][0][0]
    
    x_L = np.concatenate(Trajectory['Trajectory'][0][0]['x_L_AllDistances'])
    y_L = np.concatenate(Trajectory['Trajectory'][0][0]['y_L_AllDistances'])
    z_L = np.concatenate(Trajectory['Trajectory'][0][0]['z_L_AllDistances'])

    # yaw is a rotation around the z-axis
    T_Yaw = [[np.cos(Yaw), -np.sin(Yaw), 0],
             [np.sin(Yaw), np.cos(Yaw), 0],
             [0, 0, 1]]

    T_Pitch = [[np.cos(Yaw), 0, -np.sin(Yaw)],
               [0, 1, 0],
               [np.sin(Yaw), 0, np.cos(Yaw)]]

    T_Roll = [[1, 0, 0],
              [0, np.cos(Roll), -np.sin(Roll)],
              [0, np.sin(Roll), np.cos(Roll)]]

    T = np.dot(np.dot(T_Yaw, T_Pitch), T_Roll)

    x_R = T[0, 0]*x_L + T[0, 1]*y_L + T[0, 2]*z_L
    y_R = T[1, 0]*x_L + T[1, 1]*y_L + T[1, 2]*z_L
    z_R = T[2, 0]*x_L + T[2, 1]*y_L + T[2, 2]*z_L

    x_W = x_R + PositionLinW[0] + inputData['turbineX'][inputData['turbineLidar']]
    y_W = y_R + PositionLinW[1] + inputData['turbineY'][inputData['turbineLidar']]
    z_W = z_R + PositionLinW[2] + inputData['hubHeight'][inputData['turbineLidar']]

    for i in range(len(x_W)):
        if z_W[i] < 0:
            z_W[i] = 0.01

    return x_W, y_W, z_W

def getPointsVLOS(x_W, y_W, z_W, inputData):

    # retrieves the line of sight wind speed (vlos) from a wind field

    # input: parameter, wind field
    # output: v_los
    Parameter = scipy.io.loadmat(inputData['LidarParams'])
    Trajectory = scipy.io.loadmat(inputData['LidarTrajectory'])

    a = Parameter['Parameter'][0][0][0]['a'][0][0][0]
    nWeights = len(a)    
    nDataPoints = len(x_W)

    x_LW = np.ones(nDataPoints)*inputData['turbineX'][inputData['turbineLidar']]
    y_LW = np.ones(nDataPoints)*inputData['turbineY'][inputData['turbineLidar']]
    z_LW = np.ones(nDataPoints)*inputData['hubHeight'][inputData['turbineLidar']]

    # calculation of the normalized laser vector
    LaserVector_Wx = [np.transpose(x_W)-np.transpose(x_LW)][0]
    LaserVector_Wy = [np.transpose(y_W)-np.transpose(y_LW)][0]
    LaserVector_Wz = [np.transpose(z_W)-np.transpose(z_LW)][0]

    NormedLaserVector_Wx = np.zeros(len(LaserVector_Wx))
    NormedLaserVector_Wy = np.zeros(len(LaserVector_Wy))
    NormedLaserVector_Wz = np.zeros(len(LaserVector_Wz))

    for i in range(nDataPoints):
        NormLaserVector_W = np.sqrt( LaserVector_Wx[i]**2 + LaserVector_Wy[i]**2 + LaserVector_Wz[i]**2 )
        if NormLaserVector_W == 0:
            NormLaserVector_W = 1
        NormedLaserVector_Wx[i] = LaserVector_Wx[i]/NormLaserVector_W
        NormedLaserVector_Wy[i] = LaserVector_Wy[i]/NormLaserVector_W
        NormedLaserVector_Wz[i] = LaserVector_Wz[i]/NormLaserVector_W

    BackscatterNormedLaserVector_Wx = -NormedLaserVector_Wx
    BackscatterNormedLaserVector_Wy = -NormedLaserVector_Wy
    BackscatterNormedLaserVector_Wz = -NormedLaserVector_Wz

    # Calculation of considered points in the laser beam
    Points_WFx = np.zeros(nDataPoints*nWeights)
    Points_WFy = np.zeros(nDataPoints*nWeights)
    Points_WFz = np.zeros(nDataPoints*nWeights)
    for i in range(nDataPoints):
        #print(x_W[i]*np.ones(nWeights))
        Points_WFx[i*nWeights:((i*nWeights)+nWeights)] = x_W[i]*np.ones(nWeights) + (a*BackscatterNormedLaserVector_Wx[i]*np.ones(nWeights))
        Points_WFy[i*nWeights:((i*nWeights)+nWeights)] = y_W[i]*np.ones(nWeights) + (a*BackscatterNormedLaserVector_Wy[i]*np.ones(nWeights))
        Points_WFz[i*nWeights:((i*nWeights)+nWeights)] = z_W[i]*np.ones(nWeights) + (a*BackscatterNormedLaserVector_Wz[i]*np.ones(nWeights))

    Points_WFy = Points_WFy - inputData['turbineY'][inputData['turbineLidar']]
    Points_WFy = inputData['turbineY'][inputData['turbineLidar']] - Points_WFy

    return Points_WFx, Points_WFy, Points_WFz

def VLOS(x_W, y_W, z_W, inputData, Upts):

    Parameter = scipy.io.loadmat(inputData['LidarParams'])
    Trajectory = scipy.io.loadmat(inputData['LidarTrajectory'])

    a = Parameter['Parameter'][0][0][0]['a'][0][0][0]
    nWeights = len(a)    
    nDataPoints = len(x_W)
    f_L_d = Parameter['Parameter'][0][0][0]['f_L_d'][0][0][0]

    x_L = np.concatenate(Trajectory['Trajectory'][0][0]['x_L_AllDistances'])
    y_L = np.concatenate(Trajectory['Trajectory'][0][0]['y_L_AllDistances'])
    z_L = np.concatenate(Trajectory['Trajectory'][0][0]['z_L_AllDistances']) 
    #t = Parameter.t

    # origin of the lidar and the origin of the wind coordinate system
    x_LW = np.ones(nDataPoints)*inputData['turbineX'][inputData['turbineLidar']]
    y_LW = np.ones(nDataPoints)*inputData['turbineY'][inputData['turbineLidar']]
    z_LW = np.ones(nDataPoints)*inputData['hubHeight'][inputData['turbineLidar']]

    # calculation of the normalized laser vector (vector from the zero in the wind to the trajectory point in the wind)
    # lidar in the wind does not need to be at zero
    LaserVector_Wx = [np.transpose(x_W)-np.transpose(x_LW)][0]
    LaserVector_Wy = [np.transpose(y_W)-np.transpose(y_LW)][0]
    LaserVector_Wz = [np.transpose(z_W)-np.transpose(z_LW)][0]

    NormedLaserVector_Wx = np.zeros(len(LaserVector_Wx))
    NormedLaserVector_Wy = np.zeros(len(LaserVector_Wy))
    NormedLaserVector_Wz = np.zeros(len(LaserVector_Wz))

    for i in range(nDataPoints):
        NormLaserVector_W = np.sqrt( LaserVector_Wx[i]**2 + LaserVector_Wy[i]**2 + LaserVector_Wz[i]**2 )
        NormedLaserVector_Wx[i] = LaserVector_Wx[i]/NormLaserVector_W
        NormedLaserVector_Wy[i] = LaserVector_Wy[i]/NormLaserVector_W
        NormedLaserVector_Wz[i] = LaserVector_Wz[i]/NormLaserVector_W

    BackscatterNormedLaserVector_Wx = NormedLaserVector_Wx
    BackscatterNormedLaserVector_Wy = NormedLaserVector_Wy
    BackscatterNormedLaserVector_Wz = NormedLaserVector_Wz

    u_W = Upts
    v_W = np.zeros(nDataPoints*nWeights)
    w_W = np.zeros(nDataPoints*nWeights) 

    RelativeWindVector_W = [u_W, v_W, w_W]
    BackScatterNormed_W = [BackscatterNormedLaserVector_Wx, BackscatterNormedLaserVector_Wy, BackscatterNormedLaserVector_Wz]

    for i in range(3):
        tmp = BackScatterNormed_W[i]
        for j in range(nDataPoints):
            if j == 0:
                tmp1 = tmp[j]*np.ones(nWeights)
            else:
                tmp1 = np.concatenate((tmp1, tmp[j]*np.ones(nWeights)))  
        if i == 0:
            BackscatterNormedLaserVector_Wx = tmp1  
        elif i == 1:
            BackscatterNormedLaserVector_Wy = tmp1
        else:
            BackscatterNormedLaserVector_Wz = tmp1

    BackScatterNormed_W1 = [BackscatterNormedLaserVector_Wx, BackscatterNormedLaserVector_Wy, BackscatterNormedLaserVector_Wz]

    vlos_d = np.multiply(RelativeWindVector_W, BackScatterNormed_W1)
    v_los = np.zeros(nDataPoints)

    for i in range(nDataPoints):
        v_los[i] = np.dot(vlos_d[0, i*nWeights:((i*nWeights)+nWeights)], f_L_d)

    return v_los


def rotatedCoordinates(layout):
    # this function rotates the coordinates so that the flow direction is now
    # alligned with the x-axis. This makes computing wakes and wake overlap
    # much simpler
    windDirection = layout.windDirection
    xTurb = layout.turbineX
    yTurb = layout.turbineY
    nTurbs = len(xTurb)

    RotAng = np.radians(windDirection)
    # find the center of the xy-domain, i.e. the mean
    xCenter = np.mean(xTurb)
    yCenter = np.mean(yTurb)

#    if inputData['outputFlowField'] or inputData['visualizeHorizontal']:
#        # initialize output vectors
#        Xrot = np.zeros(X.shape)
#        Yrot = np.zeros(Y.shape)
#
#        # number of samples in the x and y coordinates
#        nSamplesX = X.shape[2]
#        nSamplesY = Y.shape[1]
#        nSamplesZ = Z.shape[0]
#
#        # demean the X and Y grids and the x and y turbine coordinates
#        # in order to rotate
#        X = X - xCenter
#        Y = Y - yCenter
#
#        # rotate the grid
#        for i in range(nSamplesZ):
#            for j in range(nSamplesY):
#                for k in range(nSamplesX):
#                    Xrot[i, j, k] = X[i, j, k]*np.cos(RotAng) - Y[i, j, k]*np.sin(RotAng) 
#                    Yrot[i, j, k] = X[i, j, k]*np.sin(RotAng) + Y[i, j, k]*np.cos(RotAng)

    xTurb = xTurb - xCenter
    yTurb = yTurb - yCenter
    xTurbRot = np.zeros(nTurbs)
    yTurbRot = np.zeros(nTurbs)

    # rotate the turbine coordinates
    for i in range(nTurbs):
        xTurbRot[i] = xTurb[i]*np.cos(RotAng) - yTurb[i]*np.sin(RotAng)
        yTurbRot[i] = xTurb[i]*np.sin(RotAng) + yTurb[i]*np.cos(RotAng)

    # put the mean back in and return the X, Y domains and the rotated turbine locations
#    if inputData['outputFlowField'] or inputData['visualizeHorizontal']:
#        return Xrot+xCenter, Yrot+yCenter, xTurbRot+xCenter, yTurbRot+yCenter
#    else:
    return xTurbRot+xCenter, yTurbRot+yCenter


def avgVelocity(X, Y, Z, Ufield, xTurb, yTurb, zTurb, D, turbI, model, cSet):

    # this function averages the velocity across the rotor.
    # The mean wind speed across the rotor is used.
    tilt = cSet.tiltAngles[turbI]
    yaw = cSet.yawAngles[turbI]
    rotorPts = int(np.round(np.sqrt(model.rotorPts)))

    xCenter = xTurb
    zR = (D/2.)*np.cos(np.radians(tilt))
    yR = (D/2.)*np.cos(np.radians(yaw))
#    xR = (D/2.)*np.sin(np.radians(yaw)) + (D/2.)*np.sin(np.radians(tilt))
    yPts = np.linspace(yTurb-yR, yTurb+yR, rotorPts)
    zPts = np.linspace(zTurb-zR, zTurb+zR, rotorPts)

    # this function averages the velocity across the whole rotor
    # number of points in the X and Y domain
    nSamplesX = X.shape[2]
    nSamplesY = Y.shape[1]
    nSamplesZ = Z.shape[0]

    # define the rotor plane along with the points associated with the rotor
    Xtmp = []
    Ytmp = []
    Ztmp = []
    Utmp = []

    # if X and Y are large, resample the domains to only include points that
    # we care about. This significantly speeds up the process
    for i in range(nSamplesZ):
        for j in range(nSamplesY):
            for k in range(nSamplesX):
                if X[i, j, k] >= (xTurb-D/4.) and X[i, j, k] <= (xTurb+D/4.):
                    dist = np.hypot(Y[i, j, k] - yTurb, Z[i, j, k] - zTurb)
                    if dist <= (D/2.):
                        Xtmp.append(X[i, j, k])
                        Ytmp.append(Y[i, j, k])
                        Ztmp.append(Z[i, j, k])
                        Utmp.append(Ufield[i, j, k])

    # interpolate the points to find the average velocity across the rotor
    if len(Xtmp) == 0 or len(Ytmp) == 0 or len(Ztmp) == 0:
        print('Too few points for outputFlowField, ' +
              'please run again with more points')

    utmp = []

    for i in range(rotorPts):
        for j in range(rotorPts):
            dist = np.hypot(yPts[i] - yTurb, zPts[j] - zTurb)
            if dist <= (D/2.):
                utmp.append(griddata((Xtmp, Ytmp, Ztmp), Utmp,
                            (xCenter, yPts[i], zPts[j]), method='nearest'))

    if model.avgCube:
        return np.mean(utmp**3)**(1./3.)
    else:
        return np.mean(utmp)


def combineWakes(Uinf, Ueff, Ufield, Uwake, model):
    # this function allows for different wake superpositions
    # freestream linear superposition
    if model.combineWakes == 0:
        vel = Uinf - ((Uinf - Uwake) + (Uinf - Ufield))

    # local velocity linear superposition
    elif model.combineWakes == 1:
        vel = Uinf - ((Ueff - Uwake) + (Uinf - Ufield))

    # sum of squares freestream superposition
    elif model.combineWakes == 2:
        vel = Uinf - np.sqrt((Uinf - Uwake)**2 + (Uinf - Ufield)**2)

    # sum of squares local velocity superposition
    elif model.combineWakes == 3:
        vel = Ueff - np.sqrt((Ueff - Uwake)**2 + (Uinf - Ufield)**2)
    else:
        raise NameError('No valid wake combination model was specified')

    return vel

def computeAddedTI(X, Y, Z, UfieldOrig, xTurb, yTurb, zTurb,
                   turbI, upwindTurbines, model, layout, cSet, output):
    # this function is necessary for computing the TI that impacts the specified turbine
    # loop through all of the turbines and compute only the velocity at the rotor of interest
    # computes areaoverlap and TI contribution
    # the TI is addition is limited to a specified distance upstream (default = 15*D)
    TI_a = model.TIa
    TI_b = model.TIb
    TI_c = model.TIc
    TI_d = model.TId

    yaw = cSet.yawAngles
    tilt = cSet.tiltAngles

    TI_0 = layout.turbulenceIntensity
    D = [turb.rotorDiameter for turb in layout.turbines]
    ke = model.wakeExpansion
    kd = model.wakeDeflection
    ad = model.ad  # lateral wake deflection a + b*X
    bd = model.bd  # lateral wake deflection a + b*X

    # determine the effective velocity generated by that specific turbine
    AreaOverlap = np.zeros(layout.nTurbs)
    TI_calc = np.zeros(layout.nTurbs)

    for turbIdx in upwindTurbines:

        Uwake = copy.copy(UfieldOrig)

        # compute the x and y offset due to yaw
        yR = (D[turbIdx]/2.)*np.cos(np.radians(yaw[turbIdx]))
        xR = (D[turbIdx]/2.)*np.sin(np.radians(yaw[turbIdx]))

        dist = np.hypot(xTurb[turbI]-xTurb[turbIdx], yTurb[turbI]-yTurb[turbIdx])
        if dist > (model.TIdistance(D[turbIdx])):
            continue

        wake = wakeModels.GAUSS(model, layout, cSet, output, turbIdx)
        # cycle through the previously defined X and Y domains
        xIdx = 0
        for yIdx in range(Y.shape[1]):
            for zIdx in range(Z.shape[0]):
                # use GAUSS wake model
                if model.WakeModel == 2:
                    Uwake[zIdx, yIdx, xIdx] = wake.V(UfieldOrig[zIdx, yIdx, xIdx], X[zIdx, yIdx, xIdx]-xTurb[turbI], Y[zIdx, yIdx, xIdx]-yTurb[turbI], Z[zIdx, yIdx, xIdx]-zTurb[turbI])

                # use FLORIS or Jensen wake model
                elif model.WakeModel == 0 or model.WakeModel == 1:    

                    if (X[zIdx, yIdx, xIdx] > (xTurb[turbIdx]-abs(xR)) and 
                        Y[zIdx, yIdx, xIdx] > (yTurb[turbIdx] - 2*D[turbIdx]) and Y[zIdx, yIdx, xIdx] < (yTurb[turbIdx] + 2*D[turbIdx]) and
                        Z[zIdx, yIdx, xIdx] > (zTurb[turbIdx] - 2*D[turbIdx]) and Z[zIdx, yIdx, xIdx] < (zTurb[turbIdx] + 2*D[turbIdx])):
                        yDisp = wakeModels.Jimenez(np.radians(yaw[turbIdx]), Ct[turbIdx], kd[turbIdx], X[zIdx, yIdx, xIdx]-xTurb[turbIdx], D[turbIdx], ad, bd)
                        zDisp = wakeModels.Jimenez(np.radians(tilt[turbIdx]), Ct[turbIdx], kd[turbIdx], X[zIdx, yIdx, xIdx]-xTurb[turbIdx], D[turbIdx], ad, bd)
                    else:
                        yDisp = 0.0
                        zDisp = 0.0
                        continue

                    # define the edges of the wake
                    xTurbY = xTurb[turbIdx] + ( (Y[zIdx, yIdx, xIdx]-yTurb[turbIdx])*np.tan(-np.radians(yaw[turbIdx])))
                    rWake = ke[turbI]*(X[zIdx, yIdx, xIdx]-(xTurb[turbI]-xR))
                    rCenterY = xTurb[turbI] + yDisp
                    rCenterZ = xTurb[turbI] + zDisp

                    # compute the velocity deficit
                    # define the edges of the wake
                    rWake = ke[turbI]*(X[zIdx, yIdx, xIdx]-(xTurb[turbI]-xR)) + (D[turbI]/2)
                    rCenterY = yTurb[turbI] + yDisp
                    rCenterZ = zTurb[turbI] + zDisp

                    rtmp = np.sqrt( (Y[zIdx, yIdx, xIdx] - rCenterY)**2 + (Z[zIdx, yIdx, xIdx] - rCenterZ)**2 )
                    if (X[zIdx, yIdx, xIdx] >= xTurb[turbI] and rtmp <= rWake):

                        # FLORIS model
                        if model.WakeModel == 1:
                            c = wakeModels.FLORIS(X[zIdx, yIdx, xIdx], Y[zIdx, yIdx, xIdx], Z[zIdx, yIdx, xIdx], yDisp, zDisp, xTurb[turbIdx], yTurb[turbIdx], zTurb[turbIdx], inputData, turbIdx)
                            velDef = Uinf*2.*a[turbIdx]*c
                        
                        # Jensen model
                        elif model.WakeModel == 0:
                            c = wakeModels.Jensen(inputData, turbIdx, X[zIdx, yIdx, xIdx], xTurb[turbIdx])
                            velDef = Uinf*2.*a[turbIdx]*c

                    else:
                        velDef = 0.0

        # compute percent overlap
        if turbIdx == turbI:
            AreaOverlap[turbIdx] = 0.0
        else:
            AreaOverlap[turbIdx] = computeOverlap(Uwake, UfieldOrig)

        # Niayifar and Porte-Agel, "A new analytical model for wind farm power prediction", 2015
        if turbI == turbIdx or xTurb[turbI] < xTurb[turbIdx]:
            TI_calc[turbIdx] = 0.0

        else:
            if (xTurb[turbI]-xTurb[turbIdx]) > 0:
                TI_calc[turbIdx] = TI_a*(output.aI[turbIdx]**TI_b)*(TI_0**TI_c)*(((xTurb[turbI]-xTurb[turbIdx])/D[turbIdx])**(TI_d))
            else:
                TI_calc[turbIdx] = 0.0
    
    
    # compute components of TI added
    TI_added = []
    TIidx = []
    Overlap = []
    for i in upwindTurbines:
        if AreaOverlap[i]*TI_calc[i] > 0.0:
            TI_added.append(AreaOverlap[i]*TI_calc[i])
            TIidx.append(i) 
            Overlap.append(AreaOverlap[i])

    return TI_added, TIidx, Overlap


def computeOverlap(Ufield, UfieldOrig):

    # compute wakeOverlap based on the number of points that are not
    # freestream velocity, i.e. affected by a wake
    idx = Ufield.shape
    numPts = idx[0]*idx[1]*idx[2]
    count = 0.
    for i in range(idx[0]):
        for j in range(idx[1]):
            for k in range(idx[2]):
                if Ufield[i, j, k] < 0.99*UfieldOrig[i, j, k]:
                    count = count + 1

    perOverlap = count/numPts
    return perOverlap

def computePower(Ueff, inputData):

    # compute the power at each turbine based on the effective velocity determined using the wake model
    nTurbs = len(inputData['turbineX'])
    rho = inputData['airDensity']
    D = inputData['rotorDiameter']
    pP = inputData['pP']
    pT = inputData['pT']
    gE = inputData['generatorEfficiency']
    CpCtTable = inputData['TurbineInfo']['CpCtWindSpeed']

    # turbine operation
    yaw = inputData['yawAngles']
    tilt = inputData['tilt']
    Ct = inputData['Ct']
    Cp = inputData['Cp']

    powerOut = np.zeros(nTurbs)
    for i in range(nTurbs):
        A = np.pi*(D[i]/2.)**2
        Cptmp = Cp[i]*(np.cos(yaw[i]*np.pi/180.)**pP[i])*(np.cos((tilt[i])*np.pi/180.)**pT[i])
        powerOut[i] = 0.5*rho*A*Cptmp*gE[i]*Ueff[i]**3

    return powerOut


def initializeFlowField(X, Y, Z, layout):

    # initialize the flow field used in the 3D model based on shear using the
    # power log law
    Ufield = np.zeros(X.shape)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            for k in range(X.shape[2]):
                Ufield[i, j, k] = layout.windSpeed * (Z[i, j, k] /
                                    layout.turbineZ[0])**layout.shear

    return Ufield

def outputUpts(inputData, X, Y, Z, Ufield):

    # compute the output velocity at the points specified in inputData
    if inputData['points']:
        xPts = inputData['xPts']
        yPts = inputData['yPts']
        zPts = inputData['zPts']
    elif inputData['Lidar']:
        x_W, y_W, z_W = getWindCoords(inputData)
        xTraj, yTraj, zTraj = getPointsVLOS(x_W, y_W, z_W, inputData)
        xPts = xTraj
        yPts = yTraj
        zPts = zTraj

    nSamplesX = X.shape[2]
    nSamplesY = Y.shape[1]
    nSamplesZ = Z.shape[0]

    Upts = np.zeros(len(xPts))
    count = 0
    for kk in range(len(xPts)):
        nextPt = 0
        for i in range(nSamplesZ):
            for j in range(nSamplesY):
                for k in range(nSamplesX):
                    if X[i, j, k] == xPts[kk] and Y[i, j, k] == yPts[kk] and Z[i, j, k] == zPts[kk]:
                        if nextPt == 0:
                            Upts[count] = Ufield[i, j, k]
                            count = count + 1
                        nextPt = 1

    return Upts
