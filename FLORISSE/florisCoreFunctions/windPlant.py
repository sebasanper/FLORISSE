# -*- coding: utf-8 -*-
import numpy as np
import copy
import outputClasses.outputs

import florisCoreFunctions.windPlantFunctions as wPFs


def windPlant(model, layout, cSet, *argv):
    # Create an output object based on the layout. This will hold the results
    # of the FLORIS run and the settings used to run it
    if argv[0]:
        output = outputClasses.outputs.fullOutput(model, layout, cSet)
    else:
        output = outputClasses.outputs.powerOutput(model, layout, cSet)

    # Rotate the frame of reference such that the wind is alligned with x-axis
    xTurb, yTurb = wPFs.rotatedCoordinates(layout)
    output.rotLocX = xTurb
    output.rotLocY = yTurb

    zTurb = np.array(layout.locZ)
    D = [turb.rotorDiameter for turb in layout.turbines]

    # Generate grid points at the swept area of every turbine
    X, Y, Z = wPFs.sweptAreaGrid(model, layout, xTurb, yTurb, zTurb)

    # Save the flowfield prediction of every turbine in Utp, the weird notation
    # is a trick for tuple concatenation
    Utp = np.zeros(X.shape + (layout.nTurbs,))

    # sort turbine coordinates from front to back
    sortedTurbIds = [i[0] for i in sorted(enumerate(xTurb), key=lambda x:x[1])]

    # initialize flow field with a uniform shear layer
    UfieldOrig = wPFs.initializeFlowField(Z, layout)
    Ufield = copy.copy(UfieldOrig)

    for turbI in sortedTurbIds:
        # compute effective wind speed at turbine by taking the average
        # velocity across the rotor disk
        output.windSpeed[turbI] = wPFs.avgVelocity(X, Y, Z, Ufield, xTurb[turbI],
                                              yTurb[turbI], zTurb[turbI],
                                              D[turbI], turbI, model, cSet)

        output = wPFs.computeCpCtPoweraI(layout, cSet, output, turbI)

        # compute the initial added turbulence at the rotor
        upWindTurbines = sortedTurbIds[:sortedTurbIds.index(turbI)]

        TI_added = wPFs.computeAddedTI(np.atleast_3d(UfieldOrig[turbI, :, :]), xTurb, yTurb, zTurb,
         Utp[turbI,:,:,:], turbI, upWindTurbines, model, layout, output)

        # add turbulence via sum of squares
        output.TI[turbI] = np.linalg.norm(TI_added + [layout.TI_0])

        # Instantiate the wake of this turbine
        output.wakes[turbI] = model.wake(model, layout, cSet, output, turbI)

        # Compute the velocity field according to this wake
        tol = (np.abs(np.sin(cSet.yawAngles[turbI])) *
               (layout.turbines[turbI].rotorDiameter)/2)
        dwDist = X[:, 0, 0]-xTurb[turbI]
        Yrel = Y - yTurb[turbI]
        Zrel = Z - zTurb[turbI]
        Utp[:, :, :, turbI] = computeVelocity(dwDist, Yrel, Zrel, UfieldOrig,
                                              output.wakes[turbI], tol)
        
        # Combine the wake of the turbine with the flowfield so far
        Ufield = model.wakeCombine(UfieldOrig, output.windSpeed[turbI],
                                   Ufield, Utp[:, :, :, turbI])

    return output


def velAtLocations(X, Y, Z, output):
    UfieldOrig = wPFs.initializeFlowField(Z, output.layout)
    Ufield = copy.copy(UfieldOrig)

    for turbI in range(output.layout.nTurbs):
        tol = (np.abs(np.sin(output.cSet.yawAngles[turbI])) *
               (output.layout.turbines[turbI].rotorDiameter)/2)
        Xrel = X[:, 0, 0]-output.rotLocX[turbI]
        Yrel = Y - output.rotLocY[turbI]
        Zrel = Z - output.layout.locZ[turbI]
        Uturb = computeVelocity(Xrel, Yrel, Zrel, UfieldOrig,
                                output.wakes[turbI], tol)
        Ufield = output.model.wakeCombine(
                UfieldOrig, output.windSpeed[turbI], Ufield, Uturb)
    return Ufield


def computeVelocity(dwDist, Y, Z, Uin, wake, tol):
    U = copy.copy(Uin)
    # cycle through the grid generated above
    for xIdx in range(dwDist.shape[0]):
        if dwDist[xIdx] >= -tol:
            yDisp, zDisp = wake.displ(dwDist[xIdx])
            yMat = Y[xIdx, :, :]-yDisp
            zMat = Z[xIdx, :, :]-zDisp
            wakeMask = wake.B(dwDist[xIdx], yMat, zMat)
            uWake = wake.V(Uin[xIdx, :, :], dwDist[xIdx], yMat, zMat)
            U[xIdx, :, :] = uWake*wakeMask + Uin[xIdx, :, :]*~wakeMask
    return U
