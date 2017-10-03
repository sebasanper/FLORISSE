# -*- coding: utf-8 -*-
import numpy as np


class callbackList:
    """A list like object that calls the supplied callback function on every
    element that is set"""
    def __init__(self, values, callback):
        self._values = values
        self.callback = callback

    def __getitem__(self, idx):
        return self._values[idx]

    def __setitem__(self, idx, value):
        self._values[idx] = value
        self.callback()

    def __len__(self):
        return len(self._values)

    def __str__(self):
        return '['+', '.join(map(str, self._values))+']'

    def __repr__(self):
        return self.__str__()


# Define a base control set with all the angles set to zero and bladepitch=1.9
class Neutral:
    """A control set for the turbines in the windfarm layout. This object holds
    yaw and tilt angles for every turbine. If required bladepitch angles for
    every turbine are also stored here.

    In addition to that the controlset holds four attributes that are derived
    from the yaw and tilt angles."""
    def __init__(self, layout):
        self.nTurbs = layout.nTurbs
        # If an element in yaw or tilt is changed the other attributes need to
        # change is well. This is accomplished by making a simple object called
        # angleList. It behaves similar as a list but calls computeAlphas
        # every time an element in the list is changed
        self._phis = [0 for i in range(self.nTurbs)]  # turbine angle [rad]

        self.wakeDir = [np.array([0, 0, 0]) for i in range(self.nTurbs)]
        self.Cvec = [np.zeros((2, 2)) for i in range(self.nTurbs)]
        self.Rvec = [np.zeros((3, 3)) for i in range(self.nTurbs)]

        self._yawAngles = callbackList([0 for i in range(self.nTurbs)],
                                       self.setDerivedAttrs)  # Yaw angle [deg]
        self.tiltAngles = [0 for i in range(self.nTurbs)]  # Tilt angle [deg]
        self.bladePitch = [1.9 if turb.usePitch else 0
                           for turb in layout.turbines]  # blade pitch [deg]

    # If either yawAngles or tiltAngles is replaced in its entirety the set
    # command is intercepted to rebind self.computeAlphas to the new list.
    @property
    def yawAngles(self):
        return self._yawAngles

    @yawAngles.setter
    def yawAngles(self, values):
        if len(values) == self.nTurbs:
            self._yawAngles = callbackList(values, self.setDerivedAttrs)
            self.setDerivedAttrs()
        else:
            raise Exception('length of yawAngles should be %d' % self.nTurbs)

    @property
    def tiltAngles(self):
        return self._tiltAngles

    @tiltAngles.setter
    def tiltAngles(self, values):
        if len(values) == self.nTurbs:
            self._tiltAngles = callbackList(values, self.setDerivedAttrs)
            self.setDerivedAttrs()
        else:
            raise Exception('length of yawAngles should be %d' % self.nTurbs)

    @property
    def phis(self):
        return self._phis

    def setDerivedAttrs(self):
        windDir = np.array([-1, 0, 0])
        for i in range(self.nTurbs):
            # Notice that tilt is defined clockwise, sign change!
            R = np.dot(self.rotMatrixZ(np.radians(self._yawAngles[i])),
                       self.rotMatrixY(np.radians(-self._tiltAngles[i])))
            thrustDir = np.dot(R, windDir)
            self._phis[i] = np.arccos(np.dot(thrustDir, windDir))
            if thrustDir[0] == -1:
                self.wakeDir[i] = thrustDir
            else:
                self.wakeDir[i] = np.array([0, thrustDir[1], thrustDir[2]])
                self.wakeDir[i] = self.wakeDir[i]/np.linalg.norm(self.wakeDir[i])
            self.Cvec[i] = np.dot(R[1:, 1:], R[1:, 1:].T)
            self.Rvec[i] = R

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


class Yawed(Neutral):
    def __init__(self, layout):
        super().__init__(layout)
        self.yawAngles[0] = 20
