# -*- coding: utf-8 -*-
from scipy.optimize import minimize, check_grad, approx_fprime
import autograd.numpy as np
from autograd import grad, value_and_grad

import florisCoreFunctions.windPlant as windPlant


class modelOptimizer(object):
    def __init__(self, model, layout, cSet, target, layoutPars = [], cSetPars = []):
        self.model, self.layout, self.cSet = model, layout, cSet
        self.target, self.layoutPars, self.cSetPars = target, layoutPars, cSetPars

        self.outputUnOptim = windPlant.windPlant(self.model, self.layout,
                                                 self.cSet, True)
        self.grad = grad(self.cost)

    def optimize(self):
        self.cost(self.x0)
        minimize(self.cost, self.x0, method='SLSQP', bounds=self.bnds,
                 jac=grad(self.cost), options={'disp': True})
        self.outputOptim = windPlant.windPlant(self.model, self.layout,
                                               self.cSet, True)
        return self.outputOptim

    def costInner(self):
        temp = []
        for keyL, valuesL in self.layoutPars.items():
            for i in range(len(valuesL)):
                for keyC, valuesC in self.cSetPars.items():
                    for ii in range(len(valuesC)):
                        setattr(self.layout, keyL, valuesL[i])
                        setattr(self.cSet, keyC, valuesC[ii])
                        output = windPlant.windPlant(self.model, self.layout, self.cSet, False)
                        temp.append(np.array(output.power) - np.array(self.target[i][ii]))
        dif = np.array(temp)
        return np.sqrt(np.sum(dif**2))


class alphaBetaOptimizer(modelOptimizer):
    def __init__(self, model, *args):
        self.x0 = [model.alpha, model.beta]
        self.bnds = [(0.0, 1.0), (0.0, 1.0)]
        super().__init__(model, *args)

    def cost(self, x):
        self.model.alpha, self.model.beta = x
        return self.costInner()


class kaKbOptimizer(modelOptimizer):
    def __init__(self, model, *args):
        self.x0 = [model.ka, model.kb, model.ad,
                   model.bd, model.aT, model.bT]
        self.bnds = [(0.0, 1.0), (np.finfo(np.float).eps, 1.0), (-10.0, 10.0), (-.5, .5), (-10.0, 10.0), (-.5, .5)]
        super().__init__(model, *args)

    def cost(self, x):
        self.model.ka, self.model.kb, self.model.ad, self.model.bd, self.model.aT, self.model.bT = x
        return self.costInner()


class TIOptimizer(modelOptimizer):
    def __init__(self, model, *args):
        self.x0 = [model.TIa, model.TIb, model.TIc, model.TId]
        self.bnds = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
        super().__init__(model, *args)

    def cost(self, x):
        self.model.TIa, self.model.TIb, self.model.TIc, self.model.TId = x
        return self.costInner()
