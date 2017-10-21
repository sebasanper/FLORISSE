# -*- coding: utf-8 -*-
import autograd.numpy as np
np.tanh.defvjp(lambda g, ans, vs, gvs, x: g * (1 - ans**2))


class Wake:
    """This class combines the information from the wake deflection and
    velocity objects to predict the position and velocity profile of a wake.
    In addition the wake boundaries at the rotor and optionally at the side can
    be defined here in a continous way using the tanh filter. This allows for
    automatic gradient generation."""

    def __init__(self, model, layout, cSet, output, turbI):
        self.wakeSlicer = model.velClass(model, layout, cSet, output, turbI).wakeSlicer
        self.displ = model.deflClass(model, layout, cSet, output, turbI).displ
        # Ri holds the inverse of rotationMatrix R
        self.Ri = cSet.Rvec[turbI].T

    def V(self, U, x, Y, Z):
        yDisp, zDisp = self.displ(x)
        wSlice = self.wakeSlicer(U, x, Y-yDisp, Z-zDisp)

        # Broadcast x to make a matrix: (x + 0*y) = np.broadcast_to(x, y.shape)
        tempCoords = np.einsum('ij, jmn -> imn',
                               self.Ri, np.stack([x + 0*Y, Y, Z]))
        rotorBoundary = self.tanHFilter(tempCoords[0], 0, 10)
        wakeBoundary = self.tanHFilter(wSlice.BT, wSlice.BV, 10)
        mask = rotorBoundary*wakeBoundary

        return wSlice.V*mask + U*(1-mask)

    def tanHFilter(self, x, loc, sharpness):
        return (1 + np.tanh(sharpness*(x-loc)))/2
