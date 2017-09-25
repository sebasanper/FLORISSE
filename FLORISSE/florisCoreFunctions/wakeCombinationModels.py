# -*- coding: utf-8 -*-
import numpy as np


# freestream linear superposition
def FLS(Uinf, Ueff, Ufield, Uwake):
    return Uinf - ((Uinf - Uwake) + (Uinf - Ufield))


# local velocity linear superposition
def LVLS(Uinf, Ueff, Ufield, Uwake):
    return Uinf - ((Ueff - Uwake) + (Uinf - Ufield))


# sum of squares freestream superposition
def SOSFS(Uinf, Ueff, Ufield, Uwake):
    return Uinf - np.sqrt((Uinf - Uwake)**2 + (Uinf - Ufield)**2)


# sum of squares local velocity superposition
def SOSLVS(Uinf, Ueff, Ufield, Uwake):
    return Ueff - np.sqrt((Ueff - Uwake)**2 + (Uinf - Ufield)**2)
