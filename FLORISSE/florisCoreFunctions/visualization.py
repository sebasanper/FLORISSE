# -*- coding: utf-8 -*-

def visualizeOutput():
        # ============================================================
    # Outputs based on input selection
    # ============================================================


    rotorPts = inputData['rotorPts']  # number of rotor points evaluated over the rotor
    rotorPts = int(np.round(np.sqrt(rotorPts)))      # evaluate the number of points in the y-z direction (equal in each direction)
    idxLidar = inputData['turbineLidar']  # turbine with the Lidar

    # Generate grid to evaluate the wake model
    nSamplesX = inputData['nSamplesX']  # number of X samples, used for visualization
    nSamplesY = inputData['nSamplesY']  # number of Y samples
    nSamplesZ = inputData['nSamplesZ']  # number of Z samples
    xLen = np.linspace(inputData['xLen'][0], inputData['xLen'][1], nSamplesX)  # x domain
    yLen = np.linspace(inputData['yLen'][0], inputData['yLen'][1], nSamplesY)  # y domain
    zLen = np.linspace(inputData['zLen'][0], inputData['zLen'][1], nSamplesZ)  # z domain
    dx = xLen[1]-xLen[0]  # grid spacing in the x direction
    dy = yLen[1]-yLen[0]  # grid spacing in the y direction

    # rotating the coordinates - flow always comes form the west
    xTurbOrig = xTurb
    yTurbOrig = yTurb
    zTurbOrig = zTurb

    # output flow field
    if inputData['outputFlowField']:
        Ueff, Ufield, xLen, yLen = velocityModel(inputData)
        outputData['Ufield'] = Ufield
        outputData['xLen'] = xLen
        outputData['yLen'] = yLen
    # output lidar parameters
    elif inputData['Lidar']:
        Ueff, Ufield, X, Y, Z, Upts, vlos = velocityModel(inputData)
        outputData['Ufield'] = Ufield
        outputData['X'] = X
        outputData['Y'] = Y
        outputData['Z'] = Z
        outputData['vlos'] = vlos
        outputData['Upts'] = Upts
    # output velocities at specified points
    elif inputData['points']:
        Ueff, Upts = velocityModel(inputData)
        outputData['Upts'] = Upts

    # different grids generated based on input selection
    # want output of the flow field or to visualize the flow
    if inputData['outputFlowField'] or inputData['visualizeHorizontal']:
        if inputData['visualizeHorizontal']:
            inputData['nSamplesZ'] = 1
            zLen = inputData['hubHeight'][0]
            
        Z, Y, X = np.meshgrid(zLen, yLen, xLen, indexing='ij')
        X, Y, xTurb, yTurb = utilities.rotatedCoordinates(inputData, X, Y, Z)

    # output generic points of the flow field - Ex: used for custom lidar model
    elif inputData['points'] or inputData['Lidar']:

        # generate compressed grid with individual x, y, z points
        X = []
        Y = []
        Z = []
        xTurb, yTurb = utilities.rotatedCoordinates(inputData, X, Y, Z)

        if inputData['points']:
            xPts = inputData['xPts']
            yPts = inputData['yPts']
            zPts = inputData['zPts']
        elif inputData['Lidar']:        
            x_W, y_W, z_W = utilities.getWindCoords(inputData)
            xTraj, yTraj, zTraj = utilities.getPointsVLOS(x_W, y_W, z_W, inputData)
            xPts = xTraj
            yPts = yTraj
            zPts = zTraj

        Xt = np.concatenate([xTurb, xPts])
        Yt = np.linspace(yTurb[0]-(D[0]/2), yTurb[0]+(D[0]/2), rotorPts)
        Zt = np.linspace(zTurb[0]-(D[0]/2), zTurb[0]+(D[0]/2), rotorPts)

        X = np.zeros((len(Zt), len(Yt), len(Xt)))
        Y = np.zeros((len(Zt), len(Yt), len(Xt)))
        Z = np.zeros((len(Zt), len(Yt), len(Xt)))

        # generate grid
        count = 0.
        for k in range(len(Xt)):
            if k < nTurbs:
                Yt = np.linspace(yTurb[k]-(D[k]/2), yTurb[k]+(D[k]/2), rotorPts)
                Zt = np.linspace(zTurb[k]-(D[k]/2), zTurb[k]+(D[k]/2), rotorPts)
            else:
                idx = k - nTurbs
                Yt = [yPts[idx]]
                Zt = [zPts[idx]]
            for j in range(len(Yt)):
                for i in range(len(Zt)):
                    X[i, j, k] = Xt[k]
                    Y[i, j, k] = Yt[j]
                    Z[i, j, k] = Zt[i]

def rotatedCoordinates(layout):
    # this function rotates the coordinates so that the flow direction is now
    # alligned with the x-axis. This makes computing wakes and wake overlap
    # much simpler
    windDirection = layout.windDirection
    xTurb = layout.locX
    yTurb = layout.locY

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
    xTurbRot = np.zeros(layout.nTurbs)
    yTurbRot = np.zeros(layout.nTurbs)

    # rotate the turbine coordinates
    for i in range(layout.nTurbs):
        xTurbRot[i] = xTurb[i]*np.cos(RotAng) - yTurb[i]*np.sin(RotAng)
        yTurbRot[i] = xTurb[i]*np.sin(RotAng) + yTurb[i]*np.cos(RotAng)

    # put the mean back in and return the X, Y domains and the rotated turbine locations
#    if inputData['outputFlowField'] or inputData['visualizeHorizontal']:
#        return Xrot+xCenter, Yrot+yCenter, xTurbRot+xCenter, yTurbRot+yCenter
#    else:
    return xTurbRot+xCenter, yTurbRot+yCenter