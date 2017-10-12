# -*- coding: utf-8 -*-
import autograd.numpy as np
np.hypot.defvjp(lambda g, ans, vs, gvs, x, y: g*x/ans)
np.hypot.defvjp(lambda g, ans, vs, gvs, x, y: g*y/ans, argnum=1)


# freestream linear superposition
def FLS(Uinf, Ueff, Ufield, Uwake):
    return Uinf - ((Uinf - Uwake) + (Uinf - Ufield))


# local velocity linear superposition
def LVLS(Uinf, Ueff, Ufield, Uwake):
    return Uinf - ((Ueff - Uwake) + (Uinf - Ufield))


# sum of squares freestream superposition
def SOSFS(Uinf, Ueff, Ufield, Uwake):
    return Uinf - np.hypot(Uinf - Uwake, Uinf - Ufield)


# sum of squares local velocity superposition
def SOSLVS(Uinf, Ueff, Ufield, Uwake):
    return Ueff - np.hypot(Ueff - Uwake, Uinf - Ufield)
