# ==============================================================================
# New analytical flow model - main file
# ==============================================================================

import numpy as np
import copy

import florisCoreFunctions.utilities as utilities
import florisCoreFunctions.wakeModels as wakeModels
import outputClasses.basic


def windPlant(model, layout, cSet):

    output = outputClasses.basic.output(layout)

    # this function uses an engineering model to compute the time
    # averaged wake of a wind turbine you can model the wake
    # using Gaussian, FLORIS, or Jensen parameters

    # wake parameters
    TI_0 = layout.turbulenceIntensity
    D = [turb.rotorDiameter for turb in layout.turbines]
    ke = model.wakeExpansion
    kd = model.wakeDeflection
    ad = model.ad  # lateral wake deflection a + b*X
    bd = model.bd  # lateral wake deflection a + b*X
    rotorPts = int(np.round(np.sqrt(model.rotorPts)))

    # turbine operation
    yaw = cSet.yawAngles
    tilt = cSet.tiltAngles

    xTurb = layout.turbineX  # x locations of the turbines
    yTurb = layout.turbineY  # y locations
    zTurb = layout.turbineZ  # z lodcations

    # generate compressed grid
    xTurb, yTurb = utilities.rotatedCoordinates(layout)

    Xt = xTurb
    Yt = np.linspace(yTurb[0]-(D[0]/2), yTurb[0]+(D[0]/2), rotorPts)
    Zt = np.linspace(zTurb[0]-(D[0]/2), zTurb[0]+(D[0]/2), rotorPts)

    X = np.zeros((len(Zt), len(Yt), len(Xt)))
    Y = np.zeros((len(Zt), len(Yt), len(Xt)))
    Z = np.zeros((len(Zt), len(Yt), len(Xt)))

    for k in range(len(Xt)):
        Yt = np.linspace(yTurb[k]-(D[k]/2), yTurb[k]+(D[k]/2), rotorPts)
        Zt = np.linspace(zTurb[k]-(D[k]/2), zTurb[k]+(D[k]/2), rotorPts)
        for j in range(len(Yt)):
            for i in range(len(Zt)):
                X[i, j, k] = Xt[k]
                Y[i, j, k] = Yt[j]
                Z[i, j, k] = Zt[i]

    # initialize z location
    zTurb = np.array(zTurb)

    # sort turbine coordinates from front to back
    sortedTurbIds = [i[0] for i in sorted(enumerate(xTurb), key=lambda x:x[1])]

    # initialize flow field with a uniform shear layer
    Ufield = utilities.initializeFlowField(X, Y, Z, layout)
    UfieldOrig = copy.copy(Ufield)

    for turbI in sortedTurbIds:
        Uwake = copy.copy(UfieldOrig)

        # compute effective wind speed at each turbine
        # take the average across the rotor disk
        output.windSpeed[turbI] = utilities.avgVelocity(X, Y, Z, Ufield,
                xTurb[turbI], yTurb[turbI], zTurb[turbI], D[turbI], turbI, model, cSet)

        # adjust Cp/Ct based on local velocity (both tables were generated using FAST, Jonkman et. al. 2005)
        # Compute Cp and Ct using wind speed and possibly pitch
        if layout.turbines[turbI].usePitch:
            output.Cp[turbI] = layout.turbines[turbI].Cp(output.windSpeed[turbI], cSet.bladePitch[turbI])
            output.Ct[turbI] = layout.turbines[turbI].Ct(output.windSpeed[turbI], cSet.bladePitch[turbI])
        else:
            output.Cp[turbI] = layout.turbines[turbI].Cp(output.windSpeed[turbI])
            output.Ct[turbI] = layout.turbines[turbI].Ct(output.windSpeed[turbI])

        if output.Ct[turbI] >= 1.0:
            output.Ct[turbI] = 0.99999999

        # Compute axial induction factor
        output.aI[turbI] = (0.5 / np.cos(yaw[turbI] * np.pi / 180.)) * (1 -
                np.sqrt(1-output.Ct[turbI]*np.cos(yaw[turbI]*np.pi/180)))

        # compute the x offset due to turbine rotation
        xR = (D[turbI]/2.)*np.sin(np.radians(yaw[turbI])) + (D[turbI]/2.)*np.sin(np.radians(tilt[turbI]))

        # compute the initial added turbulence at the rotor
        upWindTurbines = sortedTurbIds[:sortedTurbIds.index(turbI)]
        
        TI_added, TIidx, AreaOverlap = utilities.computeAddedTI(
                np.atleast_3d(X[:, :, turbI]), np.atleast_3d(Y[:, :, turbI]),
                np.atleast_3d(Z[:, :, turbI]), np.atleast_3d(UfieldOrig[:, :, turbI]),
                xTurb, yTurb, zTurb, turbI, upWindTurbines, model, layout, cSet, output)

        # add turbulence via sum of squares
        output.TI[turbI] = np.linalg.norm(TI_added + [TI_0])
        wake = wakeModels.GAUSS(model, layout, cSet, output, turbI)
#        print(dir(wake))

        # cycle through the grid generated above
        for xIdx in range(X.shape[2]):
            for yIdx in range(Y.shape[1]):
                for zIdx in range(Z.shape[0]):
                    # use Gaussian wake model
                    if model.WakeModel == 2:

                        if X[zIdx, yIdx, xIdx] >= xTurb[turbI] - D[turbI]:
                            Uwake[zIdx, yIdx, xIdx] = wake.V(UfieldOrig[zIdx, yIdx, xIdx], X[zIdx, yIdx, xIdx]-xTurb[turbI], Y[zIdx, yIdx, xIdx]-yTurb[turbI], Z[zIdx, yIdx, xIdx]-zTurb[turbI])

                    # use FLORIS or Jensen wake model
                    elif model.WakeModel == 0 or model.WakeModel == 1:    
                        if (X[zIdx, yIdx, xIdx] > (xTurb[turbI]-abs(xR)) and 
                           Y[zIdx, yIdx, xIdx] > (yTurb[turbI] - 2*D[turbI]) and 
                           Y[zIdx, yIdx, xIdx] < (yTurb[turbI] + 2*D[turbI]) and
                           Z[zIdx, yIdx, xIdx] > (zTurb[turbI] - 2*D[turbI]) and
                           Z[zIdx, yIdx, xIdx] < (zTurb[turbI] + 2*D[turbI])):
                            yDisp = wakeModels.Jimenez(np.radians(yaw[turbI]), Ct[turbI], kd[turbI], X[zIdx, yIdx, xIdx]-xTurb[turbI], D[turbI], ad, bd)
                            zDisp = wakeModels.Jimenez(np.radians(tilt[turbI]), Ct[turbI], kd[turbI], X[zIdx, yIdx, xIdx]-xTurb[turbI], D[turbI], ad, bd)
                        else:
                            yDisp = 0.0
                            zDisp = 0.0
                            continue

                        # define the edges of the wake
                        rWake = ke[turbI]*(X[zIdx, yIdx, xIdx]-(xTurb[turbI]-xR)) + (D[turbI]/2)
                        rCenterY = yTurb[turbI] + yDisp
                        rCenterZ = zTurb[turbI] + zDisp

                        rtmp = np.sqrt( (Y[zIdx, yIdx, xIdx] - rCenterY)**2 + (Z[zIdx, yIdx, xIdx] - rCenterZ)**2 )
                        if (X[zIdx, yIdx, xIdx] >= xTurb[turbI] and rtmp <= rWake):

                            # FLORIS model
                            if model.WakeModel == 1:
                                c = wakeModels.FLORIS(X[zIdx, yIdx, xIdx], Y[zIdx, yIdx, xIdx], Z[zIdx, yIdx, xIdx], yDisp, zDisp, xTurb[turbI], yTurb[turbI], zTurb[turbI], inputData, turbI)
                                velDef = layout.windSpeed*(1-2.*a[turbI]*c)

                            # Jensen model
                            elif model.WakeModel == 0:
                                c = wakeModels.Jensen(inputData, turbI, X[zIdx, yIdx, xIdx], xTurb[turbI])
                                velDef = layout.windSpeed*2.*a[turbI]*c

                        else:
                            velDef = 0.0

        Ufield = utilities.combineWakes(UfieldOrig, output.windSpeed[turbI], Ufield, Uwake, model)

    output.computePower(layout, cSet)
    return output

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
        Ueff, Ufield, xLen, yLen = WakeModel(inputData)
        outputData['Ufield'] = Ufield
        outputData['xLen'] = xLen
        outputData['yLen'] = yLen
    # output lidar parameters
    elif inputData['Lidar']:
        Ueff, Ufield, X, Y, Z, Upts, vlos = WakeModel(inputData)
        outputData['Ufield'] = Ufield
        outputData['X'] = X
        outputData['Y'] = Y
        outputData['Z'] = Z
        outputData['vlos'] = vlos
        outputData['Upts'] = Upts
    # output velocities at specified points
    elif inputData['points']:
        Ueff, Upts = WakeModel(inputData)
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




def Optimize():
    
    # =============================================================
    # perform optimization(s)
    # =============================================================
# This does not work currently
    if inputData['axial_opt']:
        CpOpt, CtOpt, bladePitchOpt = OptModules.axialOpt(inputData)
        inputData['Cp'] = CpOpt
        inputData['Ct'] = CtOpt     
        outputData['Cp_opt'] = CpOpt
        outputData['Ct_opt'] = CtOpt
        inputData['bladePitch'] = bladePitchOpt
    elif inputData['yaw_opt']:
        fileloc = []
        yawOpt = OptModules.yawOpt(inputData)
        inputData['yawAngles'] = yawOpt
        outputData['yaw_opt'] = yawOpt

    # ==============================================================
    # rerun the wake model with the optimized parameters
    # ==============================================================

    if inputData['axial_opt'] or inputData['yaw_opt']:
        Ueff = WakeModel(inputData)          

        powerOut = utilities.computePower(Ueff, inputData)
        outputData['powerOut'] = powerOut
        outputData['Ueff'] = Ueff
        outputData['yawAngles'] = inputData['yawAngles']
        outputData['bladePitch'] = inputData['bladePitch']

        powerOpt = np.sum(outputData['powerOut'])
        powerGain = 100*(powerOpt - power0)/power0
        print('Power gain = ', powerGain, '%')
        if powerGain < 0.0:
            outputData['bladePitch'] = 1.9*np.ones(len(inputData['turbineX']))
            outputData['yawAngles'] = np.zeros(len(inputData['turbineX']))
