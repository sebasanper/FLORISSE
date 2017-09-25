# -*- coding: utf-8 -*-
import numpy as np
import scipy.io


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
