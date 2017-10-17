# -*- coding: utf-8 -*-
import autograd.numpy as np
import copy
import outputClasses.outputs

import florisCoreFunctions.windPlantFunctions as wPFs
from florisCoreFunctions.wake import Wake
import pdb


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
    zTurb = layout.zLoc
    D = [turb.rotorDiameter for turb in layout.turbines]

    # Generate grid points at the swept area of every turbine
    X, Y, Z = wPFs.sweptAreaGrid(model, layout)
    # Utp will hold numpy arrays with the velocity field at each turbine
    Utp = []

    # initialize flow field with a uniform shear layer, Hub heigth of the first
    # turbine is used as the characteristic height
    UfieldOrig = wPFs.initializeFlowField(Z, layout.windSpeed, layout.shear,
                                          layout.turbines[0].hubHeight)
    Ufield = copy.copy(UfieldOrig)

    for sortedIndex in range(layout.nTurbs):
        turbI = layout.sortedTurbIds[sortedIndex]

        # compute effective wind speed at turbine by taking the average
        # velocity across the rotor disk
        output.windSpeed.append(wPFs.avgVelocity(
                X, Y, Z, Ufield, xTurb[turbI], yTurb[turbI], zTurb[turbI],
                D[turbI], cSet.yawAngles[turbI], cSet.tiltAngles[turbI], model))
        output = wPFs.computeCpCtPoweraI(layout.turbines[turbI], cSet.bladePitch[turbI], cSet.yawAngles[turbI],
                                         cSet.tiltAngles[turbI], layout.airDensity, output)

        # compute the added turbulence at the rotor cause by upwind turbines
        upWindTurbines = layout.sortedTurbIds[:sortedIndex]
        TI_added = wPFs.computeAddedTI(
                np.atleast_3d(UfieldOrig[turbI, :, :]), xTurb, yTurb, zTurb,
                Utp, turbI, upWindTurbines, model, layout, output)

        # add turbulence via sum of squares
        output.TI.append(np.linalg.norm(np.array(TI_added + [layout.TI_0])))

        # Instantiate the wake of this turbine
        output.wakes.append(Wake(model, layout, cSet, output, turbI))

        # Compute the velocity field as predicted by this wake
        tipOffset = (10 + np.sin(cSet.phis[turbI]) *
                     (layout.turbines[turbI].rotorDiameter)/2)
        dwDist = X[:, 0, 0]-xTurb[turbI]
        Yrel = Y - yTurb[turbI]
        Zrel = Z - zTurb[turbI]

        Utp.append(computeVelocity(dwDist, Yrel, Zrel, UfieldOrig,
                                   output.wakes[-1], tipOffset))

        Ufield = model.wakeCombine(UfieldOrig, output.windSpeed[sortedIndex],
                                   Ufield, Utp[sortedIndex])

    output.reorderParameters(layout.sortedTurbIds)
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
    U = []
    for xIdx in range(dwDist.shape[0]):
        if dwDist[xIdx] >= -tol:
            U.append(wake.V(Uin[xIdx, :, :], dwDist[xIdx],
                            Y[xIdx, :, :], Z[xIdx, :, :]))
        else:
            U.append(Uin[xIdx, :, :])
    return np.array(U)
