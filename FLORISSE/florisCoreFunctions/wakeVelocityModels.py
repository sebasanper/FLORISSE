# -*- coding: utf-8 -*-
import autograd.numpy as np


class WakeSlice:
    """All the wake code is vectorized in the y and z dimension. Computing
    the predicted velocity of the wake thus always happens in slices with one
    x-coordinate and a Y and Z matrix.

    These inputs than return a velocity prediction for that particular slice.
    To plot this properly and make the velocity profiles behave continous a
    wake boundary should be specified. This requires two components, a boundary
    value such as the radial distance from the centerline and a threshold.

    The wake object than combines the velocity profile and boundary information
    together with the turbine rotor the determine the wake position"""
    def __init__(self, velocity, boundaryValue, boundaryThreshold):
        self.V = velocity
        self.BV = boundaryValue
        self.BT = boundaryThreshold


class Jensen:
    """This class instantiates an object for computing the wake velocity
    profile at some point Y, Z at downwind position X accorindg to Jensen"""

    def __init__(self, model, layout, cSet, output, turbI):
        self.ke = model.wakeExpansion
        self.R = layout.turbines[turbI].rotorDiameter/2
        self.aI = output.aI[-1]

    def wakeSlicer(self, U, x, Y, Z):
        # compute the velocity based on the classic Jensen/Park model,
        # see Jensen 1983
        c = (self.R/(self.ke*(x) + self.R))**2
        v = U*(1-2.*self.aI*c)

        return WakeSlice(v, np.hypot(Y, Z), self.ke*x + self.R)


class FLORIS:
    """This class instantiates an object for computing the wake velocity
    profile at some point Y, Z at downwind position X according to the
    FLORIS model developed by gebraad et al."""

    def __init__(self, model, layout, cSet, output, turbI):
        # Save turbine specific attirbutes
        self.R = layout.turbines[turbI].rotorDiameter/2
        self.yaw = cSet.yawAngles[turbI]
        self.aI = output.aI[-1]

        # save expansion coefficient and expansion modifiers
        self.ke = model.wakeExpansion
        self.me = np.array(model.me)

        # Adjust the recovery coefficient MU for the yaw angle
        self.MU = (np.array(model.MU) /
                   (np.cos(np.radians(model.aU + model.bU*self.yaw))))

    def wakeSlicer(self, U, x, Y, Z):
        # compute the velocity based on the classic Floris model,
        # see Gebraad et al
        radius = np.hypot(Y, Z)

        # wake zone radii
        rZones = self.R + self.ke*self.me*x

        # defining wake parameters
        cZones = (self.R/(self.R + self.ke*self.MU*x))**2

        c = (((radius <= rZones[2]) ^ (radius < rZones[1]))*cZones[2] +
             ((radius <= rZones[1]) ^ (radius < rZones[0]))*cZones[1] +
             (radius <= rZones[0])*cZones[0])
        return WakeSlice(U*(1-2.*self.aI*c), radius,
                         self.R + self.ke*self.me[2]*x)


class GAUSS:
    """This class instantiates an object for computing the wake velocity
    profile at some point Y, Z at downwind position X according to the
    Porte-Age model as adapted by Jennifer Anonni"""

    def __init__(self, model, layout, cSet, output, turbI):
        self.ka = model.ka
        self.kb = model.kb
        self.alpha = model.alpha
        self.beta = model.beta

        self.veer = layout.veer
        self.D = layout.turbines[turbI].rotorDiameter
        self.Uinf = layout.windSpeed

        self.yaw = -cSet.yawAngles[turbI]  # sign reversed in literature
        self.tilt = cSet.tiltAngles[turbI]
        self.Ri = cSet.Rvec[turbI].T

        self.Ct = output.Ct[-1]
        self.TI = output.TI[-1]

    def wakeSlicer(self, U, x, Y, Z):

        # initial velocity deficits
        uR = (self.Ct*np.cos(self.yaw*np.pi/180.) /
              (2.*(1-np.sqrt(1-(self.Ct*np.cos(self.yaw*np.pi/180.))))))

        # initial Gaussian wake expansion
        sigma_z0 = self.D*0.5*np.sqrt(uR/(1 + np.sqrt(1-self.Ct)))
        sigma_y0 = (sigma_z0*(np.cos((self.yaw)*np.pi/180.)) *
                    (np.cos(self.veer*np.pi/180.)))

        # quantity that determines when the far wake starts
        x0 = (self.D*(np.cos(self.yaw*np.pi/180.) *
                      (1+np.sqrt(1-self.Ct*np.cos(self.yaw*np.pi/180.)))) /
              (np.sqrt(2)*(4*self.alpha*self.TI +
                           2*self.beta*(1-np.sqrt(1-self.Ct)))))

        # wake expansion parameters
        ky = self.ka*self.TI + self.kb
        kz = self.ka*self.TI + self.kb

        # Compute the location of the swept plane of the rotor
        xR = np.einsum('ij, jmn -> imn', self.Ri, np.stack([0*Y, Y, Z]))[0]

        sigma_y_nw = ((((x0-xR)-(x-xR))/(x0-xR))*0.501*self.D *
                      np.sqrt(self.Ct/2.) + ((x-xR)/(x0-xR))*sigma_y0)
        sigma_z_nw = ((((x0-xR)-(x-xR))/(x0-xR))*0.501*self.D *
                      np.sqrt(self.Ct/2.) + ((x-xR)/(x0-xR))*sigma_z0)
        sigma_y_fw = ky*(x - x0) + sigma_y0
        sigma_z_fw = kz*(x - x0) + sigma_z0

        # Compute the standard deviation in the near and far wake
        nwMask = x < x0
        sigma_y = sigma_y_nw*nwMask + sigma_y_fw*~nwMask
        sigma_z = sigma_z_nw*nwMask + sigma_z_fw*~nwMask

        a = ((np.cos(self.veer*np.pi/180.)**2)/(2*sigma_y**2) +
             (np.sin(self.veer*np.pi/180.)**2)/(2*sigma_z**2))
        b = (-(np.sin(2*self.veer*np.pi/180))/(4*sigma_y**2) +
             (np.sin(2*self.veer*np.pi/180.))/(4*sigma_z**2))
        c = ((np.sin(self.veer*np.pi/180.)**2)/(2*sigma_y**2) +
             (np.cos(self.veer*np.pi/180.)**2)/(2*sigma_z**2))
        totGauss = (np.exp(-(a*(Y)**2 - 2*b*(Y)*(Z) + c*(Z)**2)))

        # Compute the velocity deficit
        coreVelDef = (1-((self.Ct*np.cos(self.yaw*np.pi/180.)) /
                      (8.0*sigma_y*sigma_z/self.D**2)))
        Uwake = U-(U*(1-np.sqrt(abs(coreVelDef))))*totGauss

        return WakeSlice(Uwake, np.hypot(Y, Z), 1000)


class GAUSSThrustAngle:
    """This class instantiates an object for computing the wake velocity
    profile at some point Y, Z at downwind position X according to the
    Porte-Age model as adapted by Jennifer Anonni"""

    def __init__(self, model, layout, cSet, output, turbI):
        self.ka = model.ka
        self.kb = model.kb
        self.alpha = model.alpha
        self.beta = model.beta

        self.veer = layout.veer
        self.D = layout.turbines[turbI].rotorDiameter
        self.Ct = output.Ct[-1]
        self.yaw = cSet.yawAngles[turbI]
        self.phi = cSet.phis[turbI]
        self.wakeDir = cSet.wakeDir[turbI]
        self.C = cSet.Cvec[turbI]
        self.ellipseA = np.linalg.inv(self.C*(self.D/2)**2)

        self.C0 = 1 - np.sqrt(1 - self.Ct)
        self.E0 = self.C0**2 - 3*np.exp(1/12)*self.C0 + 3*np.exp(1/3)
        self.M0 = self.C0*(2-self.C0)
        self.sqrM0 = np.sqrt(self.M0)

        # Start of farwake
        self.x0 = (self.D*(np.cos(self.phi) *
                   (1+np.sqrt(1-self.Ct*np.cos(self.phi)))) /
                   (np.sqrt(2)*(4*model.alpha*output.TI[-1] +
                    2*model.beta*(1-np.sqrt(1-self.Ct)))))
        # Angle of near wake
        self.theta_C0 = (2*((.3*self.phi)/np.cos(self.phi)) *
                         (1-np.sqrt(1-self.Ct*np.cos(self.phi))))
        # Displacement at end of near wake
        self.deltaX0 = np.tan(self.theta_C0)*self.x0

        # Neutral covariance matrix
        self.sigNeutral_x0 = np.array([[1, 0], [0, 1]])*np.sqrt(.5)*self.D/2

        # wake expansion parameters
        self.ky = model.ka*output.TI[-1] + model.kb
        self.kz = model.ka*output.TI[-1] + model.kb

    def wakeSlicer(self, U, x, Y, Z):
        # Having yz stacked is usefull for implementing quadratic
        # multiplication of two arrays i.e. [x, y] * A * [[x, y]]
        YZ = np.stack([Y, Z])
        if x < self.x0:
            ellipse = (self.ellipseA[0][0]*Y**2 + 2*self.ellipseA[1][0]*Y*Z +
                       self.ellipseA[1][1]*Z**2)
            # Compute the standard deviation in the near and far wake
            nwCoreMask = np.sqrt(ellipse) <= (1-(x/self.x0))
            elipRatio = 1-(1-(x/self.x0))/(np.finfo(float).eps + np.sqrt(ellipse))

            S = np.linalg.inv((self.C*(self.sigNeutral_x0**2)) *
                 ((((np.finfo(float).eps + 0*(x<=0))+x*(x>0))/self.x0)**2))
            nwExp = np.exp(-0.5*np.einsum('gij, gh, hij ->ij', YZ, S, YZ) *
                           (elipRatio**2))
            nw = U*(1-self.C0*(nwCoreMask + nwExp*~nwCoreMask))
            return WakeSlice(nw, np.hypot(Y, Z), 1000)
        else:
            varWake = np.dot(self.C, (np.array([[self.ky, 0], [0, self.kz]]) *
                             (x-self.x0) + self.sigNeutral_x0)**2)
            fwExp = np.exp(-0.5*np.einsum('gij, gh, hij ->ij', YZ,
                                          np.linalg.inv(varWake), YZ))
            fw = U*(1-(1-np.sqrt(1-self.Ct *
                    np.linalg.det(np.dot(np.dot(self.C, self.sigNeutral_x0**2),
                                  np.linalg.inv(varWake)))))*fwExp)
            return WakeSlice(fw, np.hypot(Y, Z), 1000)
