# -*- coding: utf-8 -*-
import numpy as np


class jimenezDeflection:
    """This class instantiates an object for computing the downwind
    deflection of a wake according to Jimenez et al"""

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


class jimenezDeflectionThrustAngle:
    """This class instantiates an object for computing the downwind
    deflection of a wake according to Jimenez et al"""

    def __init__(self, model, layout, cSet, output, turbI):
        # Extract the model properties from model and set them in the class
        self.kd = model.wakeDeflection
        self.ad = model.ad
        self.bd = model.bd
        self.aT = model.aT
        self.bT = model.bT

        self.D = layout.turbines[turbI].rotorDiameter
        self.Ct = output.Ct[turbI]
        self.phi = cSet.phis[turbI]
        self.wakeDir = cSet.wakeDir[turbI]

        # angle of deflection
        self.xiInit = 1./2.*np.sin(self.phi)*self.Ct*np.cos(self.phi)

    def displ(self, x):
        # yaw displacement
        displ = ((self.xiInit * (15*((2*self.kd*x/self.D) + 1)**4 +
                 self.xiInit**2)/((30*self.kd/self.D)*(2*self.kd*x /
                  self.D + 1)**5.)) - (self.xiInit*self.D*(15 +
                                       self.xiInit**2.)/(30*self.kd)))
        # corrected displacement with lateral offset
        wakePos = self.wakeDir * displ

        return (wakePos[1] + (self.ad + self.bd*x),
                wakePos[2] + (self.aT + self.bT*x))


class porteAgelDeflectionThrustAngle:
    """This class instantiates an object for computing the downwind
    deflection at some downwind position X according to the
    Porte-Age model as adapted by Jennifer Anonni"""

    def __init__(self, model, layout, cSet, output, turbI):
        self.ad = model.ad
        self.bd = model.bd
        self.aT = model.aT
        self.bT = model.bT

        self.veer = layout.veer
        self.D = layout.turbines[turbI].rotorDiameter
        self.Ct = output.Ct[turbI]
        self.phi = cSet.phis[turbI]
        self.wakeDir = cSet.wakeDir[turbI]
        self.C = cSet.Cvec[turbI]

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

        self.relCoef = np.linalg.det(
                np.dot(np.dot(self.C, self.sigNeutral_x0**2), np.linalg.inv(
                        ([[self.ky, 0], [0, self.kz]]*self.M0)**2)))**0.25

    def displ(self, x):

        if x < 0:
            return 0, 0
        elif x < self.x0:
            displ = self.deltaX0*x/self.x0
        else:
            varWake = np.dot(self.C, (np.array([[self.ky, 0], [0, self.kz]]) *
                             (x-self.x0) + self.sigNeutral_x0)**2)
            lnInnerTerm = np.linalg.det(np.dot(varWake, np.linalg.inv(
                           np.dot(self.C, self.sigNeutral_x0**2))))**.25
            lnTerm = (((1.6+self.sqrM0)*(1.6*lnInnerTerm-self.sqrM0)) /
                      ((1.6-self.sqrM0)*(1.6*lnInnerTerm+self.sqrM0)))
            displ = (self.deltaX0 + (self.theta_C0*self.E0/5.2) *
                     self.relCoef*np.log(lnTerm))

        wakePos = self.wakeDir * displ

        return (wakePos[1] + (self.ad + self.bd*x),
                wakePos[2] + (self.aT + self.bT*x))
