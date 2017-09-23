# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# def visualizeHorizontal - plot the freestream flow field at hub height
# def visualizeCut - plot the freestream flow field at some x-coordinate
# def visualizeLidar - plot the output of the Stuttgart lidar model implemented
# ==============================================================================


def visualizeHorizontal(xLen, yLen, zLen, Ufield, inputData):

    # plot the freestream flow field at hub height
    print('Plotting Flow Field Horizontal...\n')

    Uinf = inputData['windSpeed']
    windDirection = inputData['windDirection']
    D = inputData['rotorDiameter']
    yaw = inputData['yawAngles']
    xTurb = inputData['turbineX']
    yTurb = inputData['turbineY']
    zTurb = inputData['turbineZ']
    turbineLidar = inputData['turbineLidar']

    # number of turbines
    nTurbs = len(xTurb)

    # rotation angle for the turbines the wind turbines
    RotAng = -np.radians(windDirection - 270)

    # plot horizontal flow field
    plt.figure(figsize=(30, 10))
    plt.contourf(xLen, yLen, Ufield[0, :, :], 50,
                 cmap='coolwarm', vmin=3.0, vmax=Uinf)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=15)

    if inputData['Lidar']:

        yLidar = inputData['yLidar']
        xLidar = inputData['xLidar']

        numLocs = len(np.unique(xLidar))
        xLocs = np.zeros(numLocs+1)
        ymin = np.zeros(numLocs+1)
        ymax = np.zeros(numLocs+1)
        cols = yLidar.columns
        ymin[0] = yTurb[turbineLidar]
        ymax[0] = yTurb[turbineLidar]
        xLocs[0] = xTurb[turbineLidar]
        for i in range(1, numLocs+1):
            xLocs[i] = xTurb[turbineLidar] + np.unique(xLidar)[i-1]
            ymin[i] = yTurb[turbineLidar] + np.min(yLidar[cols[i-1]])
            ymax[i] = yTurb[turbineLidar] + np.max(yLidar[cols[i-1]])
            plt.plot([xLocs[i], xLocs[i]], [ymin[i], ymax[i]], 'k')

        plt.plot(xLocs, ymin, 'k', linewidth=3)
        plt.plot(xLocs, ymax, 'k', linewidth=3)
        plt.plot([xTurb[turbineLidar], xTurb[turbineLidar] + 700],
                 [yTurb[turbineLidar], yTurb[turbineLidar]],
                 'k--', linewidth=3)

    for idx in range(nTurbs):

        if zTurb[idx] == zTurb[0]:

            yawR = np.radians(yaw[idx])

            # x coordinate
            xRot1 = xTurb[idx] - (-(D[idx]/2))*np.sin(RotAng+yawR)
            xRot2 = xTurb[idx] - ((D[idx]/2))*np.sin(RotAng+yawR)

            # z coordinate
            yRot1 = yTurb[idx] + (-(D[idx]/2))*np.cos(RotAng+yawR)
            yRot2 = yTurb[idx] + ((D[idx]/2))*np.cos(RotAng+yawR)

            plt.plot([xRot1, xRot2], [yRot1, yRot2], 'k', linewidth=3)

    plt.axis('equal')
    plt.xlabel('x(m)', fontsize=15)
    plt.ylabel('y(m)', fontsize=15)
    plt.tick_params(which='both', labelsize=15)
    plt.title('Horizontal', fontsize=15)


def visualizeCut(xLen, yLen, zLen, Ufield, inputData):

    # plot the freestream flow field at some x-coordinate
    print('Plotting Flow Field Cut Through Slice at ',
          inputData['downLocs'], '...\n')

    Uinf = inputData['windSpeed'][0]
    windDirection = inputData['windDirection'][0]
    tilt = -1 * inputData['tilt']
    yaw = inputData['yawAngles']
    D = inputData['rotorDiameter']
    HH = inputData['hubHeight']

    xTurb = inputData['turbineX']
    yTurb = inputData['turbineY']
    zTurb = inputData['turbineZ']

    # number of turbines
    nTurbs = len(xTurb)

    # rotation angle for the turbines the wind turbines
    RotAng = -(windDirection - 270) * np.pi/180

    # plot horizontal flow field
    plt.figure(figsize=(10, 7))
    plt.contourf(-yLen, zLen, Ufield[:, :, 0], 50,
                 cmap='coolwarm', vmin=4.0, vmax=Uinf)
    plt.colorbar()

    # plot rotor
    yCirc = np.linspace(-D[0]/2, D[0]/2, 100)
    zCirc1 = np.sqrt((D[0]/2)**2 - (yCirc-yTurb[0])**2) + zTurb[0]
    zCirc2 = -np.sqrt((D[0]/2)**2 - (yCirc-yTurb[0])**2) + zTurb[0]
    plt.plot(yCirc, zCirc1, 'k', linewidth=3)
    plt.plot(yCirc, zCirc2, 'k', linewidth=3)

    # for i in range(len(xLen)):
    #   for j in range(len(yLen)):
    #       plt.plot(xLen[i],yLen[j],'k.',linewidth=3)

    # plt.axis('equal')
    plt.xlabel('y(m)', fontsize=15)
    plt.ylabel('z(m)', fontsize=15)
    strTitle = 'Cut Through Slice at ', inputData['downLocs'][0]/D[0], 'D'
    plt.title(strTitle, fontsize=15)
    plt.xlim([-2*D[0], 2*D[0]])
    plt.ylim([0.0, 2*HH[0]])
    plt.axis('equal')


def visualizeLidar(xTurb, yTurb, X, Y, Z, Ufield, inputData, vlos):

    # plot the output of the Stuttgart lidar model implemented
    Uinf = inputData['windSpeed']

    # parameters from inputData
    idxTurb = inputData['turbineLidar']
    zTurb = inputData['turbineZ'][idxTurb]

    # Lidar parameters
    xLocs = np.unique(inputData['xLidar']) + xTurb
    numLocs = len(xLocs)

    # plot Lidar panels
    plt.figure(figsize=(20, 3))
    cols = inputData['yLidar'].columns

    for kk in range(numLocs):
        plt.subplot(1, numLocs, kk+1)
        yPlot = np.unique(inputData['yLidar'][cols[kk]]) + yTurb
        zPlot = np.unique(inputData['zLidar'][cols[kk]]) + zTurb

        Uplot = np.zeros((len(yPlot), len(zPlot)))

        # find the Ufield points that correspond to yPlot and zPlot
        for ii in range(len(zPlot)):
            for jj in range(len(yPlot)):
                if (X[ii, jj, kk] == xLocs[kk] and
                        Y[ii, jj, kk] == yPlot[jj] and
                        Z[ii, jj, kk] == zPlot[ii]):
                    Uplot[ii, jj] = Ufield[ii, jj, kk]

        plt.contourf(-yPlot, zPlot, Uplot, 50,
                     cmap='coolwarm', vmin=2.0, vmax=Uinf)
        plt.xlim([-60 - yTurb, 60 - yTurb])
        plt.ylim([40, 160])
        plt.colorbar()
