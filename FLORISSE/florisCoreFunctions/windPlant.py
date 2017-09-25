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

        TI_added = wPFs.computeAddedTI(np.atleast_3d(UfieldOrig[:, :, turbI]), xTurb, yTurb, zTurb,
         Utp[:,:,turbI,:], turbI, upWindTurbines, model, layout, output)

        # add turbulence via sum of squares
        output.TI[turbI] = np.linalg.norm(TI_added + [layout.turbulenceIntensity])

        # Instantiate the wake of this turbine
        output.wakes[turbI] = model.wake(model, layout, cSet, output, turbI)

        # cycle through the grid generated above
        for xIdx in range(X.shape[2]):
            yDisp, zDisp = output.wakes[turbI].displ(X[0, 0, xIdx])
            for yIdx in range(Y.shape[1]):
                for zIdx in range(Z.shape[0]):
                    # Update the gridpoint behind the turbine and inside the wake
                    if (X[zIdx, yIdx, xIdx] >= xTurb[turbI] - D[turbI] and
                            output.wakes[turbI].B(X[zIdx, yIdx, xIdx]-xTurb[turbI], Y[zIdx, yIdx, xIdx]-yTurb[turbI]-yDisp, Z[zIdx, yIdx, xIdx]-zTurb[turbI]-zDisp)):
                        Utp[zIdx, yIdx, xIdx, turbI] = output.wakes[turbI].V(UfieldOrig[zIdx, yIdx, xIdx], X[zIdx, yIdx, xIdx]-xTurb[turbI], Y[zIdx, yIdx, xIdx]-yTurb[turbI]-yDisp, Z[zIdx, yIdx, xIdx]-zTurb[turbI]-zDisp)

        Ufield = model.wakeCombine(UfieldOrig, output.windSpeed[turbI],
                                   Ufield, Utp[:, :, :, turbI])

    return output


def velAtLocations(X, Y, Z, output):
    xTurb = output.rotLocX
    yTurb = output.rotLocY
    zTurb = output.layout.locZ
    D = [turb.rotorDiameter for turb in output.layout.turbines]

    UfieldOrig = wPFs.initializeFlowField(Z, output.layout)
    Ufield = copy.copy(UfieldOrig)

    for turbI in range(output.layout.nTurbs):
        UfieldOrigEmpty = copy.copy(UfieldOrig)
        # cycle through the grid generated above
        for xIdx in range(X.shape[2]):
            if X[0, 0, xIdx] >= xTurb[turbI]- D[turbI]:
                yDisp, zDisp = output.wakes[turbI].displ(X[0, 0, xIdx])
                dwDist = X[0, 0, xIdx]-xTurb[turbI]
                yMat = Y[:, :, xIdx]-yTurb[turbI]-yDisp
                zMat = Z[:, :, xIdx]-zTurb[turbI]-zDisp
                wakeMask = output.wakes[turbI].B(dwDist, yMat, zMat)
                uWake = output.wakes[turbI].V(UfieldOrig[:, :, xIdx], dwDist, yMat, zMat)

                UfieldOrigEmpty[:, :, xIdx] = uWake*wakeMask + UfieldOrigEmpty[:, :, xIdx]*~wakeMask

        Ufield = output.model.wakeCombine(UfieldOrig, output.windSpeed[turbI],
                                           Ufield, UfieldOrigEmpty)
    return Ufield
