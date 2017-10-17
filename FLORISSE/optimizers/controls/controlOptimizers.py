# optimization modules

from scipy.optimize import minimize
from autograd import grad

import florisCoreFunctions.windPlant as windPlant


def axialCost(x, strOpt, model, layout, cSet):
    cSet.bladePitch = x
    output = windPlant.windPlant(model, layout, cSet, False)
    return -1*sum(output.power)


def yawCost(x, strOpt, model, layout, cSet):
    cSet.yawAngles = x
    output = windPlant.windPlant(model, layout, cSet, False)
    return -1*sum(output.power)


def axialOpt(model, layout, cSet):
    x0 = cSet.bladePitch
    print('================================================================')
    print('Optimizing axial induction control...')
    print('Number of parameters to optimize = ', len(x0))
    print('================================================================\n')

    outputUnOptim = windPlant.windPlant(model, layout, cSet, True)
    powerOld = sum(outputUnOptim.power)

    bnds = [turb.betaLims for turb in layout.turbines]
    resPlant = minimize(axialCost, x0, args=('axial', model, layout, cSet),
                        method='SLSQP', bounds=bnds, options={'ftol': 0.01, 'eps': 0.5})

    cSet.bladePitch = resPlant.x
    outputOpt = windPlant.windPlant(model, layout, cSet, True)

    print('Optimal pitch angles for:')
    for i in range(layout.nTurbs):
        print('Turbine ', i, ' pitch angle = ', resPlant.x[i])

    powerGain = 100*(sum(outputOpt.power) - powerOld)/powerOld
    print('Power gain = ', powerGain, '%\n')

    return outputOpt


def yawOpt(model, layout, cSet):
    x0 = cSet.yawAngles
    print('================================================================')
    print('Optimizing wake redirection control...')
    print('Number of parameters to optimize = ', len(x0))
    print('================================================================\n')

    # put bounds on design variables
    minYaw = [-25]
    maxYaw = [25.0]

    bnds = []
    for i in range(layout.nTurbs):
        if len(minYaw) > 1:
            bnds.append((minYaw[i], maxYaw[i]))
        else:
            bnds.append((minYaw, maxYaw))
    outputUnOptim = windPlant.windPlant(model, layout, cSet, True)
    powerOld = sum(outputUnOptim.power)

    resPlant = minimize(yawCost, x0, args=('yaw', model, layout, cSet),
                        method='SLSQP', bounds=bnds, jac=grad(yawCost),
                        options={'ftol': 0.1, 'eps': 5.0})

    cSet.yawAngles = resPlant.x
    outputOpt = windPlant.windPlant(model, layout, cSet, True)

    print('Optimal yaw angles for:')
    for i in range(layout.nTurbs):
        print('Turbine ', i, ' yaw angle = ', resPlant.x[i])

    powerGain = 100*(sum(outputOpt.power) - powerOld)/powerOld
    print('Power gain = ', powerGain, '%\n')

    return outputOpt
