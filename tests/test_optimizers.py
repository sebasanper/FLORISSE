# -*- coding: utf-8 -*-
import unittest
from autograd.util import nd, check_equivalent
from autograd import value_and_grad, grad

import florisCoreFunctions.windPlant as windPlant
import inputClasses.layouts as layouts
import inputClasses.controlSettings
import inputClasses.modelData


def check_grads(fun, *args, **kwargs):
    if not args:
        raise Exception("No args given")

    exact = tuple([grad(fun, i)(*args) for i in range(len(args))])
    args = [float(x) if isinstance(x, int) else x for x in args]
    numeric = nd(fun, *args)
    check_equivalent(exact, numeric,
                     kwargs.get("rtol", 1e-4), kwargs.get("atol", 1e-6))


class TestGradients(unittest.TestCase):

    def setUp(self):
        # Select a velocity, deflection and wake summing model
        self.model = inputClasses.modelData.modelData(2, 1, 2)
        # Select wind farm layout and specify controlset
        self.layout = layouts.Layout2by2(True)
        self.layout.windDirection = -25
        self.cSet = inputClasses.controlSettings.Neutral(self.layout)

    """Tests for autograd gradient wrt yaw."""
    def test_gradient_wrt_yaw(self):
        def yawCost(x):
            self.cSet.yawAngles = x
            output = windPlant.windPlant(self.model, self.layout, self.cSet, False)
            return -1*sum(output.power)

        check_grads(yawCost, [0.0, 0.0, 0.0, 0.0], atol=1e-5)
        check_grads(yawCost, [10.0, 0.1, 3.4, 0.0])

if __name__ == '__main__':
    unittest.main()
