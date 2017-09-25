# -*- coding: utf-8 -*-

def Optimize():
    
    # =============================================================
    # perform optimization(s)
    # =============================================================
# This does not work currently
    if inputData['axial_opt']:
        CpOpt, CtOpt, bladePitchOpt = OptModules.axialOpt(inputData)
        inputData['Cp'] = CpOpt
        inputData['Ct'] = CtOpt     
        outputData['Cp_opt'] = CpOpt
        outputData['Ct_opt'] = CtOpt
        inputData['bladePitch'] = bladePitchOpt
    elif inputData['yaw_opt']:
        fileloc = []
        yawOpt = OptModules.yawOpt(inputData)
        inputData['yawAngles'] = yawOpt
        outputData['yaw_opt'] = yawOpt

    # ==============================================================
    # rerun the wake model with the optimized parameters
    # ==============================================================

    if inputData['axial_opt'] or inputData['yaw_opt']:
        Ueff = velocityModel(inputData)          

        powerOut = utilities.computePower(Ueff, inputData)
        outputData['powerOut'] = powerOut
        outputData['Ueff'] = Ueff
        outputData['yawAngles'] = inputData['yawAngles']
        outputData['bladePitch'] = inputData['bladePitch']

        powerOpt = np.sum(outputData['powerOut'])
        powerGain = 100*(powerOpt - power0)/power0
        print('Power gain = ', powerGain, '%')
        if powerGain < 0.0:
            outputData['bladePitch'] = 1.9*np.ones(len(inputData['turbineX']))
            outputData['yawAngles'] = np.zeros(len(inputData['turbineX']))
