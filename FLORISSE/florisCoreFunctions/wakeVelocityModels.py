# -*- coding: utf-8 -*-
import numpy as np


class Jensen:
    """This class instantiates an object for computing the wake velocity
    profile at some point Y, Z at downwind position X accorindg to Jensen"""

    def __init__(self, model, layout, cSet, output, turbI):
        self.ke = model.wakeExpansion
        self.R = layout.turbines[turbI].rotorDiameter/2
        self.aI = output.aI[turbI]

    def V(self, U, x, y, z):
        # compute the velocity based on the classic Jensen/Park model,
        # see Jensen 1983
        c = (self.R/(self.ke*(x) + self.R))**2
        v = U*(1-2.*self.aI*c)
        return v

    def B(self, x, y, z):
        return (np.hypot(y, z) < (self.ke*x + self.R)) & (x > 0)


class FLORIS:
    """This class instantiates an object for computing the wake velocity
    profile at some point Y, Z at downwind position X according to the
    FLORIS model developed by gebraad et al."""

    def __init__(self, model, layout, cSet, output, turbI):
        self.ke = model.wakeExpansion
        self.me = np.array(model.me)

        self.R = layout.turbines[turbI].rotorDiameter/2
        self.aI = output.aI[turbI]
        self.yaw = cSet.yawAngles[turbI]

        # Adjust the expansion coefficient MU for the yaw angle
        self.MU = (np.array(model.MU) /
                   (np.cos(np.radians(model.aU + model.bU*self.yaw))))

    def V(self, U, x, y, z):
        # compute the velocity based on the classic Floris model,
        # see Gebraad et al
        r = np.hypot(y, z)

        # wake zone radii
        rZones = self.R + self.ke*self.me*x

        # defining wake parameters
        cZones = (self.R/(self.R + self.ke*self.MU*x))**2

        c = (((r <= rZones[2]) ^ (r < rZones[1]))*cZones[2] +
             ((r <= rZones[1]) ^ (r < rZones[0]))*cZones[1] +
             (r <= rZones[0])*cZones[0])
        return U*(1-2.*self.aI*c)

    def B(self, x, y, z):
        return (np.hypot(y, z) < (self.R + self.ke*self.me[2]*x)) & (x > 0)


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
        self.aI = output.aI[turbI]
        self.Ct = output.Ct[turbI]
        self.TI = output.TI[turbI]
        self.yaw = -cSet.yawAngles[turbI]  # sign reversed in literature
        self.tilt = cSet.tiltAngles[turbI]

    def V(self, U, x, y, z):

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
        xR = y*np.tan(self.yaw*np.pi/180.)
        behindSwept = x > xR

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
        totGauss = (np.exp(-(a*(y)**2 - 2*b*(y)*(z) + c*(z)**2)))

        # Compute the velocity deficit
        coreVelDef = (1-((self.Ct*np.cos(self.yaw*np.pi/180.)) /
                      (8.0*sigma_y*sigma_z/self.D**2)))
        Uwake = U-(U*(1-np.sqrt(abs(coreVelDef))))*totGauss

        return U*~behindSwept + Uwake*behindSwept

    def B(self, x, y, z):
        return np.ones(y.shape, dtype=bool)


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
        self.Ct = output.Ct[turbI]
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
                   (np.sqrt(2)*(4*model.alpha*output.TI[turbI] +
                    2*model.beta*(1-np.sqrt(1-self.Ct)))))
        # Angle of near wake
        self.theta_C0 = (2*((.3*self.phi)/np.cos(self.phi)) *
                         (1-np.sqrt(1-self.Ct*np.cos(self.phi))))
        # Displacement at end of near wake
        self.deltaX0 = np.tan(self.theta_C0)*self.x0

        # Neutral covariance matrix
        self.sigNeutral_x0 = np.array([[1, 0], [0, 1]])*np.sqrt(.5)*self.D/2

        # wake expansion parameters
        self.ky = model.ka*output.TI[turbI] + model.kb
        self.kz = model.ka*output.TI[turbI] + model.kb

    def V(self, U, x, y, z):
        xR = y*np.tan(-self.yaw*np.pi/180.)
        if x < self.x0 and (x > xR).any():
            ellipse = (self.ellipseA[0][0]*y**2 + 2*self.ellipseA[1][0]*y*z +
                       self.ellipseA[1][1]*z**2)
            # Compute the standard deviation in the near and far wake
            nwCoreMask = np.sqrt(ellipse) <= (1-(x/self.x0)) 
            elipRatio = 1-(1-(x/self.x0))/(np.finfo(float).eps + np.sqrt(ellipse))
    
            S = np.linalg.inv((self.C*(self.sigNeutral_x0**2)) *
                 ((((np.finfo(float).eps + 0*(x<=0))+x*(x>0))/self.x0)**2))
            nwExp = (np.exp(-0.5*(elipRatio**2)*np.squeeze(np.matmul(np.matmul(
                     np.expand_dims(np.stack((y, z), 2), 2), S),
                     np.expand_dims(np.stack((y, z), 2), 3)),(2, 3))))
            nw = U*(1-self.C0*(nwCoreMask + nwExp*~nwCoreMask)*(x > xR))
            return nw
        elif x >= self.x0:
            varWake = np.dot(self.C, (np.array([[self.ky, 0], [0, self.kz]]) *
                          (x-self.x0) + self.sigNeutral_x0)**2)
            fwExp = (np.exp(-0.5*np.squeeze(np.matmul(np.matmul(
                     np.expand_dims(np.stack((y, z), 2), 2), np.linalg.inv(varWake)),
                     np.expand_dims(np.stack((y, z), 2), 3)),(2, 3))))
            fw = U*(1-(1-np.sqrt(1-self.Ct*np.cos(self.phi) *
                  np.linalg.det(np.dot(np.dot(self.C, self.sigNeutral_x0**2),
                  np.linalg.inv(varWake)))))*fwExp)
            return fw
        else:
            return U

    def B(self, x, y, z):
        return np.ones(y.shape, dtype=bool)
