# Wake models

import numpy as np


class Jimenez:
    """This class instantiates an object for computing the downwind
    deflection of a wake according to Jensen et al"""

    def __init__(self, yaw, Ct, D, model):
        # Extract the model properties from model and set them in the class
        self.kd = model.kd
        self.ad = model.ad
        self.bd = model.bd
        self.yaw = yaw
        self.Ct = Ct
        self.D = D

        # this function defines the angle at which the wake deflects in
        # relation to the yaw of the turbine this is coded as defined in the
        # Jimenez et. al. paper

        # angle of deflection
        self.xi_init = 1./2.*np.cos(self.yaw)*np.sin(self.yaw)*self.Ct
        # xi = the angle at location x, this expression is not further used,
        # yYaw uses the second order taylor expansion of xi.
        # xi = (xi_init)/(( 1 + 2*kd*(x/D) )**2)

    def displ(self, x):
        # yaw displacement
        yYaw = ((self.xi_init * (15*((2*self.kd*x/self.D) + 1)**4 +
                     self.xi_init**2)/((30*self.kd/self.D)*(2*self.kd*x /
                     self.D + 1)**5.)) - (self.xi_init*self.D*(15 +
                     self.xi_init**2.) / (30*self.kd)))

        # corrected yaw displacement with lateral offset
        return yYaw + (self.ad + self.bd*x)


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
        uR = (U*self.Ct*np.cos(self.yaw*np.pi/180.) /
              (2.*(1-np.sqrt(1-(self.Ct*np.cos(self.yaw*np.pi/180.))))))
        u0 = U*np.sqrt(1-self.Ct)

        # initial Gaussian wake expansion
        sigma_z0 = self.D*0.5*np.sqrt(uR/(U + u0))
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

        C0 = 1 - u0/self.Uinf
        M0 = C0*(2-C0)
        E0 = C0**2 - 3*np.exp(1./12.)*C0 + 3*np.exp(1./3.)

        # yaw parameters (skew angle and distance from centerline)
        theta_c0 = (2*((.3*self.yaw*np.pi/180.)/np.cos(self.yaw*np.pi/180.)) *
                    (1-np.sqrt(1-self.Ct*np.cos(self.yaw*np.pi/180.))))
        theta_z0 = (2*((.3*self.tilt*np.pi/180.)/np.cos(self.tilt*np.pi/180.))*
                    (1-np.sqrt(1-self.Ct*np.cos(self.tilt*np.pi/180.))))
        delta0 = np.tan(theta_c0)*(x0)
        delta_z0 = np.tan(theta_z0)*(x0)

        # COMPUTE VELOCITY DEFICIT
        xR = y*np.tan(self.yaw*np.pi/180.)

        if x > xR:
            if x >= (x0):
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
                          sigma_z0/(ky*kz*M0))*np.log(ln_deltaNum/ln_deltaDen)+
                          (self.aT + self.bT*x))
            else:
                sigma_y = ((((x0-xR)-(x-xR))/(x0-xR))*0.501*self.D *
                           np.sqrt(self.Ct/2.) + ((x-xR)/(x0-xR))*sigma_y0)
                sigma_z = ((((x0-xR)-(x-xR))/(x0-xR))*0.501*self.D *
                           np.sqrt(self.Ct/2.) + ((x-xR)/(x0-xR))*sigma_z0)
                delta = ((x-xR)/(x0-xR))*delta0 + (self.ad + self.bd*x)
                deltaZ = ((x-xR)/(x0-xR))*delta_z0 + (self.aT + self.bT*x)

            a = ((np.cos(self.veer*np.pi/180.)**2)/(2*sigma_y**2) +
                 (np.sin(self.veer*np.pi/180.)**2)/(2*sigma_z**2))
            b = (-(np.sin(2*self.veer*np.pi/180))/(4*sigma_y**2) +
                 (np.sin(2*self.veer*np.pi/180.))/(4*sigma_z**2))
            c = ((np.sin(self.veer*np.pi/180.)**2)/(2*sigma_y**2) +
                 (np.cos(self.veer*np.pi/180.)**2)/(2*sigma_z**2))
            totGauss = (np.exp(-(a*(y-delta)**2 - 2*b*(y-delta)*(z-deltaZ) +
                        c*(z-deltaZ)**2)))

            Unew = (U-(U*(1-np.sqrt(1-((self.Ct*np.cos(self.yaw*np.pi/180.)) /
                    (8.0*sigma_y*sigma_z/self.D**2)))))*totGauss)
        else:
            Unew = U

        return Unew
