



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
        return np.hypot(y, z) < (self.ke*x + self.R)


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
        return np.hypot(y, z) < (self.R + self.ke*self.me[2]*x)


class GAUSS:
    """This class instantiates an object for computing the wake velocity
    profile at some point Y, Z at downwind position X according to the
    Porte-Age model as adapted by Jennifer Anonni"""

    def __init__(self, model, layout, cSet, output, turbI):
        self.ka = model.ka
        self.kb = model.kb
        self.alpha = model.alpha
        self.beta = model.beta
        self.ad = model.ad
        self.bd = model.bd
        self.aT = model.aT
        self.bT = model.bT

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

        # COMPUTE VELOCITY DEFICIT
        xR = y*np.tan(self.yaw*np.pi/180.)

        behindSwept = x > xR
#        if x>xR:
        sigma_y_nw = ((((x0-xR)-(x-xR))/(x0-xR))*0.501*self.D *
                      np.sqrt(self.Ct/2.) + ((x-xR)/(x0-xR))*sigma_y0)
        sigma_z_nw = ((((x0-xR)-(x-xR))/(x0-xR))*0.501*self.D *
                      np.sqrt(self.Ct/2.) + ((x-xR)/(x0-xR))*sigma_z0)
        sigma_y_fw = ky*(x - x0) + sigma_y0
        sigma_z_fw = kz*(x - x0) + sigma_z0

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

        Uwake = (U-(U*(1-np.sqrt(1-((self.Ct*np.cos(self.yaw*np.pi/180.)) /
                 (8.0*sigma_y*sigma_z/self.D**2)))))*totGauss)
        Uwake[np.isnan(Uwake)] = 0
#        else:
#            Uwake = U
#        return Uwake
        
        return U*~behindSwept + Uwake*behindSwept

    def B(self, x, y, z):
        return np.ones(y.shape, dtype=bool)
