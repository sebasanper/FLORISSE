# -*- coding: utf-8 -*-
import numpy as np


class jimenezDeflection:
    """This class instantiates an object for computing the downwind
    deflection of a wake according to Jensen et al"""

    def __init__(self, model, layout, cSet, output, turbI):
        # Extract the model properties from model and set them in the class
        self.kd = model.wakeDeflection
        self.ad = model.ad
        self.bd = model.bd
        self.D = layout.turbines[turbI].rotorDiameter
        self.Ct = output.Ct[turbI]
        self.yaw = cSet.yawAngles[turbI]  # sign reversed in literature
        self.tilt = cSet.tiltAngles[turbI]

        # angle of deflection
        self.xiInitYaw = 1./2.*np.cos(self.yaw)*np.sin(self.yaw)*self.Ct
        self.xiInitTilt = 1./2.*np.cos(self.tilt)*np.sin(self.tilt)*self.Ct
        # xi = the angle at location x, this expression is not further used,
        # yYaw uses the second order taylor expansion of xi.
        # xiYaw = (xiInitYaw)/(( 1 + 2*kd*(x/D) )**2)

    def displ(self, x):
        # yaw displacement
        displYaw = ((self.xiInitYaw * (15*((2*self.kd*x/self.D) + 1)**4 +
                    self.xiInitYaw**2)/((30*self.kd/self.D)*(2*self.kd*x /
                     self.D + 1)**5.)) - (self.xiInitYaw*self.D*(15 +
                                          self.xiInitYaw**2.)/(30*self.kd)))
        # corrected yaw displacement with lateral offset
        displYawTotal = displYaw + (self.ad + self.bd*x)

        displTilt = ((self.xiInitTilt * (15*((2*self.kd*x/self.D) + 1)**4 +
                     self.xiInitTilt**2)/((30*self.kd/self.D)*(2*self.kd*x /
                      self.D + 1)**5.)) - (self.xiInitTilt*self.D*(15 +
                                           self.xiInitTilt**2.)/(30*self.kd)))
        displTiltTotal = displTilt + (self.ad + self.bd*x)

        return displYawTotal, displTiltTotal


class porteAgelDeflection:
    """This class instantiates an object for computing the downwind
    deflection at some downwind position X according to the
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

    def displ(self, x):

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

        C0 = 1 - np.sqrt(1 - self.Ct)
        M0 = C0*(2-C0)
        E0 = C0**2 - 3*np.exp(1./12.)*C0 + 3*np.exp(1./3.)

        # yaw parameters (skew angle and distance from centerline)
        theta_c0 = (2*((.3*self.yaw*np.pi/180)/np.cos(self.yaw*np.pi/180)) *
                    (1-np.sqrt(1-self.Ct*np.cos(self.yaw*np.pi/180.))))
        theta_z0 = (2*((.3*self.tilt*np.pi/180)/np.cos(self.tilt*np.pi/180)) *
                    (1-np.sqrt(1-self.Ct*np.cos(self.tilt*np.pi/180.))))
        delta0 = np.tan(theta_c0)*(x0)
        delta_z0 = np.tan(theta_z0)*(x0)

        if x < 0:
            delta = 0
            deltaZ = 0
        elif x < x0:
            delta = (x/x0)*delta0 + (self.ad + self.bd*x)
            deltaZ = (x/x0)*delta_z0 + (self.aT + self.bT*x)
        else:
            sigma_y = ky*(x - x0) + sigma_y0
            sigma_z = kz*(x - x0) + sigma_z0
            ln_deltaNum = ((1.6+np.sqrt(M0))*(1.6*np.sqrt(sigma_y*sigma_z /
                           (sigma_y0*sigma_z0)) - np.sqrt(M0)))
            ln_deltaDen = ((1.6-np.sqrt(M0))*(1.6*np.sqrt(sigma_y*sigma_z /
                           (sigma_y0*sigma_z0)) + np.sqrt(M0)))
            delta = (delta0 + (theta_c0*E0/5.2)*np.sqrt(sigma_y0 *
                     sigma_z0/(ky*kz*M0))*np.log(ln_deltaNum/ln_deltaDen) +
                     (self.ad + self.bd*x))
            deltaZ = (delta_z0 + (theta_z0*E0/5.2)*np.sqrt(sigma_y0 *
                      sigma_z0/(ky*kz*M0))*np.log(ln_deltaNum/ln_deltaDen) +
                      (self.aT + self.bT*x))

        return delta, deltaZ
