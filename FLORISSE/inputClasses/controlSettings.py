# -*- coding: utf-8 -*-
import numpy as np


# Define a base control set with all the angles set to zero and bladepitch=1.9
class neutral:
    """A control set for the turbines in the windfarm layout"""
    def __init__(self, layout):
        self.nTurbs = layout.nTurbs
        zeroList = [0 for i in range(self.nTurbs)]
        # If an element in yaw or tilt is changed the other attributes need to
        # change is well. This is accomplished by making a simple object called
        # angleList. It behaves similar as a list but calls computeAlphas
        # every time an element in the list is changed
        self._yawAngles = angleList(zeroList, self.computeAlphas)  # [deg]
        self._tiltAngles = angleList(zeroList, self.computeAlphas)  # [deg]
        self._alphas = zeroList  # [rad]

        self.wakeDir = [np.array([0, 0, 0]) for i in range(self.nTurbs)]
        self.Cvec = [np.array([[0, 0], [0, 0]]) for i in range(self.nTurbs)]
        self.bladePitch = [1.9 if turb.usePitch else 0
                           for turb in layout.turbines]  # [deg]

    # If either yawAngles or tiltAngles is replaced in its entirety the set
    # command is intercepted to rebind self.computeAlphas.
    @property
    def yawAngles(self):
        return self._yawAngles

    @yawAngles.setter
    def yawAngles(self, value):
        self._yawAngles = angleList(value, self.computeAlphas)
        self.computeAlphas()

    @property
    def tiltAngles(self):
        return self._tiltAngles

    @tiltAngles.setter
    def tiltAngles(self, value):
        self._tiltAngles = angleList(value, self.computeAlphas)
        self.computeAlphas()

    @property
    def alpha(self):
        return self._alphas

    def computeAlphas(self):
        print('tumtumtum')
        windDir = np.array([-1, 0, 0])
        for i in range(self.nTurbs):
            R = np.dot(self.rotMatrixZ(np.radians(self._yawAngles[i])),
                       self.rotMatrixY(np.radians(self._tiltAngles[i])))
            thrustDir = np.dot(R, windDir)
            self._alphas[i] = np.arccos(np.dot(thrustDir, windDir))
            if thrustDir[0] == -1:
                self.wakeDir = thrustDir
            else:
                self.wakeDir = np.array([0, thrustDir[1], thrustDir[2]])
                self.wakeDir = self.wakeDir/np.linalg.norm(self.wakeDir)

    # Define a rotation matrix around Z
    def rotMatrixZ(self, gamma):
        return np.array([[np.cos(gamma), -np.sin(gamma), 0],
                        [np.sin(gamma),  np.cos(gamma), 0],
                        [0, 0, 1]])

    # Define a rotation matrix around Y
    def rotMatrixY(self, tau):
        return np.array([[np.cos(tau), 0, np.sin(tau)],
                        [0, 1, 0],
                        [-np.sin(tau), 0, np.cos(tau)]])


class yawed(neutral):
    def __init__(self, layout):
        super().__init__(layout)
        self.yawAngles[0] = 20


class angleList:
    """A list like object that calls the supplied callback function on every
    element that is set"""
    def __init__(self, angles, callback):
        self._x = angles
        self.callback = callback

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x = value
        self.callback()

    def __getitem__(self, idx):
        return self._x[idx]

    def __setitem__(self, idx, value):
        self._x[idx] = value
        self.callback()
