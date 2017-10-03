# -*- coding: utf-8 -*-
import numpy as np
import copy
import outputClasses.outputs

import florisCoreFunctions.windPlantFunctions as wPFs
from florisCoreFunctions.wake import Wake


def windPlant(model, layout, cSet, *argv):
    # Create an output object based on the layout. This will hold the results
    # of the FLORIS run and the settings used to run it
    if argv[0]:
        # fullOutput saves the model layout and controlset as attributes
        output = outputClasses.outputs.fullOutput(
                  copy.copy(model), copy.copy(layout), copy.copy(cSet))
    else:
        output = outputClasses.outputs.powerOutput(model, layout, cSet)

    # Rotate the frame of reference such that the wind is alligned with x-axis
    xTurb = layout.xLocRot
    yTurb = layout.yLocRot
    zTurb = np.array(layout.zLoc)
    D = [turb.rotorDiameter for turb in layout.turbines]

    # Generate grid points at the swept area of every turbine
    X, Y, Z = wPFs.sweptAreaGrid(model, layout)

    # Save the flowfield prediction of every turbine in Utp, the weird notation
    # is a trick for tuple concatenation
    Utp = np.zeros(X.shape + (layout.nTurbs,))

    # sort turbine coordinates from front to back
    sortedTurbIds = [i[0] for i in sorted(enumerate(xTurb), key=lambda x:x[1])]

    # initialize flow field with a uniform shear layer, Hub heigth of the first
    # turbine is used as the characteristic height
    UfieldOrig = wPFs.initializeFlowField(Z, layout.windSpeed ,layout.shear,
                                          layout.turbines[0].hubHeight)
    Ufield = copy.copy(UfieldOrig)

    for turbI in sortedTurbIds:
        # compute effective wind speed at turbine by taking the average
        # velocity across the rotor disk
        output.windSpeed[turbI] = wPFs.avgVelocity(
                X, Y, Z, Ufield, xTurb[turbI], yTurb[turbI], zTurb[turbI],
                D[turbI], turbI, model, cSet)

        output = wPFs.computeCpCtPoweraI(layout, cSet, output, turbI)

        # compute the added turbulence at the rotor cause by upwind turbines
        upWindTurbines = sortedTurbIds[:sortedTurbIds.index(turbI)]
        TI_added = wPFs.computeAddedTI(
                np.atleast_3d(UfieldOrig[turbI, :, :]), xTurb, yTurb, zTurb,
                Utp[turbI, :, :, :], turbI, upWindTurbines, model, layout, output)

        # add turbulence via sum of squares
        output.TI[turbI] = np.linalg.norm(TI_added + [layout.TI_0])

        # Instantiate the wake of this turbine
        output.wakes[turbI] = Wake(model, layout, cSet, output, turbI)

        # Compute the velocity field as predicted by this wake
        tipOffset = (10 + np.sin(cSet.phis[turbI]) *
                     (layout.turbines[turbI].rotorDiameter)/2)
        dwDist = X[:, 0, 0]-xTurb[turbI]
        Yrel = Y - yTurb[turbI]
        Zrel = Z - zTurb[turbI]
        Utp[:, :, :, turbI] = computeVelocity(dwDist, Yrel, Zrel, UfieldOrig,
                                              output.wakes[turbI], tipOffset)

        # Combine the wake of the turbine with the flowfield so far
        Ufield = model.wakeCombine(UfieldOrig, output.windSpeed[turbI],
                                   Ufield, Utp[:, :, :, turbI])

    return output


def velAtLocations(X, Y, Z, output):
    UfieldOrig = wPFs.initializeFlowField(Z, output.layout.windSpeed,
                                          output.layout.shear,
                                          output.layout.turbines[0].hubHeight)
    Ufield = copy.copy(UfieldOrig)

    for turbI in range(output.layout.nTurbs):
        tipOffset = (10 + np.sin(output.cSet.phis[turbI]) *
                     (output.layout.turbines[turbI].rotorDiameter)/2)
        Xrel = X[:, 0, 0]-output.layout.xLocRot[turbI]
        Yrel = Y - output.layout.yLocRot[turbI]
        Zrel = Z - output.layout.zLoc[turbI]
        Uturb = computeVelocity(Xrel, Yrel, Zrel, UfieldOrig,
                                output.wakes[turbI], tipOffset)
        Ufield = output.model.wakeCombine(
                UfieldOrig, output.windSpeed[turbI], Ufield, Uturb)
    return Ufield


def computeVelocity(dwDist, Y, Z, Uin, wake, tol):
    U = copy.copy(Uin)
    # cycle through the grid generated above
    for xIdx in range(dwDist.shape[0]):
        if dwDist[xIdx] >= -tol:
            U[xIdx, :, :] = wake.V(Uin[xIdx, :, :], dwDist[xIdx],
                                   Y[xIdx, :, :], Z[xIdx, :, :])
    return U
