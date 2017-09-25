# -*- coding: utf-8 -*-
import numpy as np
from scipy.interpolate import griddata

# 1. rotatedCoordinates - Align the x-axis with the wind direction.
# 2. sweptAreaGrid - Generate grid points for swept areas of all the turbines
# 3. initializeFlowField - initialize the flow field used in the 3D model based
#                          on shear using the power log law
# 4. avgVelocity - compute the average velocity across a rotor
# 5. computeCpCtPoweraI - compute the power, aI and thrust + power coefficient
# 6. computeAddedTI - compute the turbulence intensity added to the wake based
#                     on turbine operation and ambient turbulence intensity
# 7. computeOverlap - compute overlap of an upstream turbine wake on swept area


def rotatedCoordinates(layout):
    # this function rotates the coordinates so that the flow direction is now
    # alligned with the x-axis. This makes computing wakes and wake overlap
    # much simpler
    windDirection = layout.windDirection
    xTurb = layout.locX
    yTurb = layout.locY

    RotAng = np.radians(windDirection)

    xTurbRot = np.zeros(layout.nTurbs)
    yTurbRot = np.zeros(layout.nTurbs)

    # rotate the turbine coordinates
    for i in range(layout.nTurbs):
        xTurbRot[i] = xTurb[i]*np.cos(RotAng) - yTurb[i]*np.sin(RotAng)
        yTurbRot[i] = xTurb[i]*np.sin(RotAng) + yTurb[i]*np.cos(RotAng)

    # Return the X and Y domains after adding the mean back to the locations
    return xTurbRot-min(xTurbRot), yTurbRot-min(yTurbRot)


def sweptAreaGrid(model, layout, xTurb, yTurb, zTurb):
    D = [turb.rotorDiameter for turb in layout.turbines]
    rotorPts = int(np.round(np.sqrt(model.rotorPts)))

    X = np.zeros([len(xTurb), rotorPts, rotorPts])
    Y = np.zeros(X.shape)
    Z = np.zeros(X.shape)

    Xt = xTurb
    for i in range(len(Xt)):
        Yt = np.linspace(yTurb[i]-(D[i]/2), yTurb[i]+(D[i]/2), rotorPts)
        Zt = np.linspace(zTurb[i]-(D[i]/2), zTurb[i]+(D[i]/2), rotorPts)
        for j in range(len(Yt)):
            for k in range(len(Zt)):
                X[i, j, k] = Xt[i]
                Y[i, j, k] = Yt[j]
                Z[i, j, k] = Zt[k]

    return X, Y, Z


def initializeFlowField(Z, layout):
    # initialize the flow field used in the 3D model based on shear using the
    # power log law
    Ufield = (layout.windSpeed * (Z/np.mean(layout.locZ))**layout.shear)

    return Ufield


def avgVelocity(X, Y, Z, Ufield, xTurb, yTurb, zTurb, D, turbI, model, cSet):
    # this function averages the velocity across the rotor.
    # The mean wind speed across the rotor is used.
    tilt = cSet.tiltAngles[turbI]
    yaw = cSet.yawAngles[turbI]
    rotorPts = int(np.round(np.sqrt(model.rotorPts)))

    zR = (D/2.)*np.cos(np.radians(tilt))
    yR = (D/2.)*np.cos(np.radians(yaw))

    yPts = np.linspace(yTurb-yR, yTurb+yR, rotorPts)
    zPts = np.linspace(zTurb-zR, zTurb+zR, rotorPts)

    # define the rotor plane along with the points associated with the rotor
    Xtmp = []
    Ytmp = []
    Ztmp = []
    Utmp = []

    for i in range(X.shape[0]):
        for j in range(Y.shape[1]):
            for k in range(Z.shape[2]):
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
                            (xTurb, yPts[i], zPts[j]), method='nearest'))

    if model.avgCube:
        return np.mean(utmp**3)**(1./3.)
    else:
        return np.mean(utmp)


def computeCpCtPoweraI(layout, cSet, output, turbI):
    # adjust Cp/Ct based on local velocity (both tables were generated
    # using FAST, Jonkman et. al. 2005)
    # Compute Cp and Ct using wind speed and possibly pitch
    turbine = layout.turbines[turbI]
    yaw = cSet.yawAngles[turbI]
    tilt = cSet.tiltAngles[turbI]

    if layout.turbines[turbI].usePitch:
        output.Cp[turbI] = turbine.Cp(output.windSpeed[turbI],
                                      cSet.bladePitch[turbI])
        output.Ct[turbI] = turbine.Ct(output.windSpeed[turbI],
                                      cSet.bladePitch[turbI])
    else:
        output.Cp[turbI] = turbine.Cp(output.windSpeed[turbI])
        output.Ct[turbI] = turbine.Ct(output.windSpeed[turbI])

    if output.Ct[turbI] >= 1.0:
        output.Ct[turbI] = 0.99999999

    Cptmp = (output.Cp[turbI]*(np.cos(yaw*np.pi/180.)**turbine.pP) *
             (np.cos(tilt*np.pi/180.)**turbine.pT))
    output.power[turbI] = (.5*layout.airDensity *
                           (np.pi*(turbine.rotorDiameter/2)**2)*Cptmp *
                           turbine.gE*output.windSpeed[turbI]**3)
    # Compute axial induction factor
    output.aI[turbI] = ((0.5 / np.cos(cSet.yawAngles[turbI] * np.pi / 180.)) *
     (1 - np.sqrt(1-output.Ct[turbI]*np.cos(cSet.yawAngles[turbI]*np.pi/180))))
    return output


def computeAddedTI(UfieldOrig, xTurb, yTurb, zTurb, Utp,
                   turbI, upwindTurbines, model, layout, output):
    # this function computes TI contribution from every upwind turbine to the
    # turbine 'turbI'
    # upwind turbine contribution is limited to a specified distance upstream
    # (default = 15*D)

    D = [turb.rotorDiameter for turb in layout.turbines]

    # determine the effective velocity generated by that specific turbine
    AreaOverlap = np.zeros(layout.nTurbs)
    TI_calc = np.zeros(layout.nTurbs)

    for turbIdx in upwindTurbines:
        if (np.hypot(xTurb[turbI]-xTurb[turbIdx],
            yTurb[turbI]-yTurb[turbIdx])>(model.TIdistance(D[turbIdx]))):
            continue

        # compute percent overlap
        Uwake = np.atleast_3d(Utp[:, :, turbIdx])
        AreaOverlap[turbIdx] = computeOverlap(Uwake, UfieldOrig)

        # Niayifar and Porte-Agel, "A new analytical model for wind farm power
        # prediction", 2015
        # TODO: discuss with Jen, originial uses TI_0 instead of TI_turb
        if (xTurb[turbI]-xTurb[turbIdx]) > 0:
            TI_calc[turbIdx] = (model.TIa*(output.aI[turbIdx]**model.TIb) *
                                (layout.TI_0**model.TIc) *
                                (((xTurb[turbI]-xTurb[turbIdx]) /
                                 D[turbIdx])**(model.TId)))
        else:
            TI_calc[turbIdx] = 0.0

    # compute components of TI added
    TI_added = []
    for i in upwindTurbines:
        if AreaOverlap[i]*TI_calc[i] > 0.0:
            TI_added.append(AreaOverlap[i]*TI_calc[i])

    return TI_added


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
