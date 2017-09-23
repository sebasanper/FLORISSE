# optimization modules

import numpy as np

from scipy.optimize import minimize

import florisCoreFunctions.main as main

def optPlant(x,strOpt,inputData):

	nTurbs 	= len(inputData['turbineX'])
	D 		= inputData['rotorDiameter'][0]
	Uinf 	= inputData['windSpeed']
	Ueff 	= main.WakeModel(inputData)

	if strOpt == 'axial':
		inputData['bladePitch'] = x
	elif strOpt == 'yaw':
		if inputData['WakeModel'] == 2:
			# due to sign convention
			inputData['yawAngles'] = x
		else:
			inputData['yawAngles'] = x

	Ueff 		= main.WakeModel(inputData)
	powerOut 	= utilities.computePower(Ueff,inputData)

	# maximize the power output of the wind farm
	powerOpt 	= -1*sum(powerOut)

	return powerOpt

def axialOpt(inputData):

	Uinf = inputData['windSpeed']
	nTurbs = len(inputData['turbineX'])

	x0 = inputData['bladePitch'] 

	# to do axial induction control, you need pre-generated cp/ct tables based on pitch and wind speed
	if inputData['TurbineInfo']['PitchCpCt'] == True:
		fCp,fCt,beta 	= inputData['TurbineInfo']['CpCtPitch']
	else:
		print('Cannot perform axial induction control without Cp/Ct tables based on pitch')
		return inputData['Cp'], inputData['Ct'], inputData['bladePitch']

	print('=====================================================================')
	print('Optimizing axial induction control...')
	print('Number of parameters to optimize = ', len(x0))
	print('=====================================================================')

	# put bounds on design variables
	bnds = []
	minbeta = np.min(beta)
	maxbeta = np.max(beta)
	for i in range(nTurbs):
		bnds.append((minbeta,maxbeta))
	
	resPlant = minimize(optPlant,x0,args=('axial',inputData),method='SLSQP',bounds=bnds,options={'ftol':0.01,'eps':0.5})

	CpOpt = []
	CtOpt = []
	bladePitchOpt = resPlant.x

	for i in range(nTurbs):
		CpOpt.append(fCp(Uinf,resPlant.x[i]))
		CtOpt.append(fCt(Uinf,resPlant.x[i]))

	return CpOpt,CtOpt,bladePitchOpt

def yawOpt(inputData):

	yaw 	= inputData['yawAngles']
	minYaw 	= [0.0]
	maxYaw 	= [25.0] 

	nTurbs 	= len(inputData['turbineX'])

	x0 		= inputData['yawAngles'] 

	print('=====================================================================')
	print('Optimizing wake redirection control...')
	print('Number of parameters to optimize = ', len(x0))
	print('=====================================================================')

	# put bounds on design variables
	bnds = []
	for i in range(nTurbs):
		if len(minYaw) > 1:
			bnds.append((minYaw[i],maxYaw[i]))
		else:
			bnds.append((minYaw,maxYaw))
	
	resPlant = minimize(optPlant,x0,args=('yaw',inputData),method='SLSQP',bounds=bnds,options={'ftol':0.1,'eps':5.0})

	yawOpt = resPlant.x

	print('Optimal yaw angles for:')
	for i in range(nTurbs):
		print('Turbine ', i, ' yaw angle = ', resPlant.x[i])

	return yawOpt
