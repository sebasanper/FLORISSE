# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# visualizeHorizontal - plot the freestream flow field at hub height
# visualizeCut - plot the freestream flow field at some x-coordinate
# visualizeLidar - plot the output of the Stuttgart lidar model implemented
# ==============================================================================


def visualizeHorizontal(flowData, output):
    # Make some of the used variables less verbose
    xTurb = output.layout.xLocRot
    yTurb = output.layout.yLocRot
    zTurb = output.layout.zLoc
    D = [turb.rotorDiameter for turb in output.layout.turbines]

    # plot horizontal flow field
    plt.figure(figsize=(30, 10))
    plt.contourf(flowData.X[:, :, 0], flowData.Y[:, :, 0], flowData.U[:, :, 0],
                 50, cmap='coolwarm', vmin=3.0, vmax=output.layout.windSpeed)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=15)

    for idx in range(output.layout.nTurbs):
        if zTurb[idx] == zTurb[0]:
            yawR = np.radians(output.cSet.yawAngles[idx])

            # x coordinate
            xRot1 = xTurb[idx] - (-(D[idx]/2))*np.sin(yawR)
            xRot2 = xTurb[idx] - ((D[idx]/2))*np.sin(yawR)

            # z coordinate
            yRot1 = yTurb[idx] + (-(D[idx]/2))*np.cos(yawR)
            yRot2 = yTurb[idx] + ((D[idx]/2))*np.cos(yawR)

            plt.plot([xRot1, xRot2], [yRot1, yRot2], 'k', linewidth=3)

    plt.axis('equal')
    plt.xlabel('x(m)', fontsize=15)
    plt.ylabel('y(m)', fontsize=15)
    plt.tick_params(which='both', labelsize=15)
    plt.title('Horizontal', fontsize=15)


def visualizeCut(flowData, output, x, turbI):
    # Make some of the used variables less verbose
    yTurb = output.layout.yLocRot
    zTurb = output.layout.zLoc
    D = [turb.rotorDiameter for turb in output.layout.turbines]

    # plot horizontal flow field
    plt.figure(figsize=(10, 7))
    plt.contourf(flowData.Y[0, :, :], flowData.Z[0, :, :], flowData.U[0, :, :],
                 50, cmap='coolwarm', vmin=4.0, vmax=output.layout.windSpeed)
    plt.colorbar()

    # plot rotor
    yCirc = np.linspace(-D[turbI]/2, D[turbI]/2, 100)
    zCirc1 = np.sqrt((D[turbI]/2)**2 - (yCirc-yTurb[turbI])**2) + zTurb[turbI]
    zCirc2 = -np.sqrt((D[turbI]/2)**2 - (yCirc-yTurb[turbI])**2) + zTurb[turbI]
    plt.plot(yCirc, zCirc1, 'k', linewidth=3)
    plt.plot(yCirc, zCirc2, 'k', linewidth=3)

    plt.xlabel('y(m)', fontsize=15)
    plt.ylabel('z(m)', fontsize=15)
    strTitle = 'Cut Through Slice at ', x/D[turbI], 'D'
    plt.title(strTitle, fontsize=15)

    plt.xlim([-2*D[turbI], 2*D[turbI]])
    plt.ylim([0.0, 2*zTurb[turbI]])
#    print(dir(plt.Axes))
#    plt.Axes.set_aspect(aspect=1)  # TODO fix aspectratio


def visualizeLidar(flowData, output, turbI):
    plt.figure(figsize=(20, 3))
    for i in range(flowData.X.shape[0]):
        plt.subplot(1, 5, i+1)
        strLoc = 'Downwind distance = %02d' % flowData.X[i, 0, 0]

        plt.contourf(flowData.Y[i, :, :], flowData.Z[i, :, :], flowData.U[i, :, :],
                     50, cmap='coolwarm', vmin=2.0, vmax=output.layout.windSpeed)
        plt.title(strLoc, fontsize=12)
        plt.xlabel('y (m)', fontsize=12)
        if i == 0:
            plt.ylabel('z (m)', fontsize=12)
        plt.xlim([-100, 100])
        plt.ylim([30, 150])
        plt.colorbar()
