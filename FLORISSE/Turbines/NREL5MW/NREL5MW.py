# default NREL 5MW turbine file
# note similar files can be produced for other turbines
import os
import numpy as np
from pandas import read_csv
from scipy.interpolate import interp1d, interp2d


class NREL5MWTurbine:
    """A wind turbine definition"""
    def __init__(self, usePitch, useTSR):
        self.rotorDiameter = 126.0
        self.hubHeight = 90
        self.NumBlades = 3
        self.pP = 1.88
        self.pT = 2.07
        self.generatorEfficiency = 1.0
        self.eta = 0.768

        self.useTSR = useTSR
        self.usePitch = usePitch
        if usePitch:
            self.Cp, self.Ct = CpCtpitch()
        else:
            self.Cp, self.Ct = CpCtWs()


def CpCtWs():
    CP = [0.        , 0.15643578, 0.31287155, 0.41306749, 0.44895632,
          0.46155227, 0.46330747, 0.46316077, 0.46316077, 0.46280642,
          0.45223111, 0.39353012, 0.3424487 , 0.2979978 , 0.25931677,
          0.22565665, 0.19636572, 0.17087684, 0.1486965 , 0.12939524,
          0.11259934, 0.0979836 , 0.08526502, 0.07419736, 0.06456631,
          0.05618541, 0.04889237, 0.]

    CT = [1.10610965, 1.09515807, 1.0227122 , 0.9196487 , 0.85190470,
          0.80328229, 0.76675469, 0.76209299, 0.76209299, 0.75083241,
          0.67210674, 0.52188504, 0.43178758, 0.36443258, 0.31049874,
          0.26696686, 0.22986909, 0.19961578, 0.17286245, 0.15081457,
          0.13146666, 0.11475968, 0.10129584, 0.0880188 , 0.07746819,
          0.06878621, 0.05977061, 0.]

    wind_speed = [0.        ,  2.5       ,  3.52338654,  4.57015961,
                  5.61693268,  6.66370575,  7.71047882,  8.75725189,
                  9.80402496, 10.85079803, 11.70448774, 12.25970155,
                  12.84125247, 13.45038983, 14.08842222, 14.75672029,
                  15.45671974, 16.18992434, 16.95790922, 17.76232421,
                  18.60489742, 19.48743891, 20.41184461, 21.38010041,
                  22.39428636, 23.45658122, 24.56926707, 30.]

    fCpInterp = interp1d(wind_speed, CP)
    fCtInterp = interp1d(wind_speed, CT)

    def fCp(Ws):
        return np.array(max(CP)) if Ws < min(wind_speed) else fCpInterp(Ws)

    def fCt(Ws):
        return np.array(.99) if Ws < min(wind_speed) else fCtInterp(Ws)

    return fCp, fCt


def CpCtpitch():
    # lookup table for Cp/Ct based on FAST

    dir_path = os.path.dirname(os.path.realpath(__file__))
    filenameCp = dir_path + '/cpPitch.csv'
    filenameCt = dir_path + '/ctPitch.csv'

    cp = read_csv(filenameCp)
    ct = read_csv(filenameCt)

    cp = cp.rename(columns={'Unnamed: 0': 'beta'})
    ct = ct.rename(columns={'Unnamed: 0': 'beta'})

    # interpolation functions for Cp and Ct
    windSpeeds = []
    beta = []
    for i in range(len(cp.columns)-1):
        windSpeeds.append(float(cp.columns[i+1]))
    for i in range(len(cp['beta'])):
        beta.append(cp['beta'][i])

    Cp = []
    Ct = []
    for i in range(len(windSpeeds)):
        for j in range(len(beta)):
            Cp.append(cp[cp.columns[i+1]][j])
            Ct.append(ct[ct.columns[i+1]][j])

    fCp = interp2d(windSpeeds, beta, Cp, kind='cubic')
    fCt = interp2d(windSpeeds, beta, Ct, kind='cubic')

    return fCp, fCt
