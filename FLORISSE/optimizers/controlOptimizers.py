# -*- coding: utf-8 -*-
from scipy.optimize import minimize
from autograd import grad

import florisCoreFunctions.windPlant as windPlant


class controlOptimizer:
    def __init__(self, model, layout, cSet):
        self.model, self.layout, self.cSet = model, layout, cSet
        self.outputUnOptim = windPlant.windPlant(self.model, self.layout,
                                                 self.cSet, True)
        self.grad = grad(self.cost)

    def optimize(self):
        minimize(self.cost, self.x0, method='SLSQP', bounds=self.bnds,
                 jac=self.grad, options={'ftol': 0.01, 'eps': 0.5})
        self.outputOptim = windPlant.windPlant(self.model, self.layout,
                                               self.cSet, True)
        powerGain = 100*(sum(self.outputOptim.power) /
                         sum(self.outputUnOptim.power) - 1)
        print('Power gain = ', powerGain, '%\n')
        return self.outputOptim


class axialOptimizer(controlOptimizer):
    def __init__(self, model, layout, cSet):
        self.x0 = cSet.bladePitch
        self.bnds = [turb.betaLims for turb in layout.turbines]
        super().__init__(model, layout, cSet)

    def cost(self, x):
        self.cSet.bladePitch = x
        output = windPlant.windPlant(self.model, self.layout, self.cSet, False)
        return -1*sum(output.power)


class yawOptimizer(controlOptimizer):
    def __init__(self, model, layout, cSet):
        self.x0 = cSet.yawAngles
        self.bnds = [(-25.0, 25.0) for i in range(layout.nTurbs)]
        super().__init__(model, layout, cSet)

    def cost(self, x):
        self.cSet.yawAngles = x
        output = windPlant.windPlant(self.model, self.layout, self.cSet, False)
        return -1*sum(output.power)


class tiltOptimizer(controlOptimizer):
    def __init__(self, model, layout, cSet):
        self.x0 = cSet.tiltAngles
        self.bnds = [(-25.0, 25.0) for i in range(layout.nTurbs)]
        super().__init__(model, layout, cSet)

    def cost(self, x):
        self.cSet.tiltAngles = x
        output = windPlant.windPlant(self.model, self.layout, self.cSet, False)
        return -1*sum(output.power)


class yawAndAxialOptimizer(controlOptimizer):
    def __init__(self, model, layout, cSet):
        self.x0 = [a for a in cSet.yawAngles] + cSet.bladePitch
        self.bnds = ([(-25.0, 25.0) for i in range(layout.nTurbs)] +
                     [turb.betaLims for turb in layout.turbines])
        super().__init__(model, layout, cSet)

    def cost(self, x):
        self.cSet.yawAngles = x[:self.layout.nTurbs]
        self.cSet.bladePitch = x[self.layout.nTurbs:]
        output = windPlant.windPlant(self.model, self.layout, self.cSet, False)
        return -1*sum(output.power)
