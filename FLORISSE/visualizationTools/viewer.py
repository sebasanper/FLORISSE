# -*- coding: utf-8 -*-
import numpy as np
import sys
from pandas import read_csv

import florisCoreFunctions.windPlant
import visualizationTools.MPLbased.MPLVisualizations as MPLViews
import visualizationTools.VTKbased.cutoffVisualizationTool as cutoffTool
import visualizationTools.VTKbased.slicerTool as slicerTool
from visualizationTools.VTKbased.writeVtiFlowfield import vtiFlowfieldWriter
from PyQt5.QtWidgets import QApplication


class flowData:
    def __init__(self, X, Y, Z, U):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.U = U


class viewer:
    def __init__(self, output):
        self.xSlice = []
        self.turbSlice = []
        self.prepared = []
        self.vtiFile = []
        self.output = output

    def showView(self, visTool, *argv):
        if visTool not in self.prepared:
            print('preparing the flowfield for plotting...\n')
            if visTool == 0:
                self.hubField = self.generateHubFlowfield()
                self.prepared.append(0)
            if visTool == 1:
                try:  # TODO: Make this catch more clear
                    self.sliceField = self.generateSliceFlowfield(*argv)
                except:
                    print('This function plots a y-z slice x meters behing' +
                          'turbine I, please specify a valid x and I')
                    return
                self.prepared.append(1)  # TODO: unset this for new settings
            if visTool == 2:
                try:  # TODO: Make this catch more clear
                    self.lidarField = self.generateLidarFlowfield(*argv)
                except:
                    print('This function plots a lidar field behind' +
                          'turbine I, please specify a valid I')
                    return
                self.prepared.append(1)  # TODO: unset this at new turbine
            if ((visTool == 3 or visTool == 4) and not
                    (self.vtiFile and self.output.prPath)):
                self.vtiFile = (self.output.prPath +
                                '/visualizationData/currentFlowData.vti')
                vFW = vtiFlowfieldWriter(self.vtiFile)
                vFW.write(self.generateFullFlowfield())
                self.prepared.extend([3, 4])

        if visTool == 0:
            MPLViews.visualizeHorizontal(self.hubField, self.output)
        elif visTool == 1:
            MPLViews.visualizeCut(self.sliceField, self.output, *argv)
        elif visTool == 2:
            MPLViews.visualizeLidar(self.lidarField, self.output, *argv)
        elif (visTool == 3 or visTool == 4) and (visTool in self.prepared):
            self.QTapp = QApplication.instance()
            if self.QTapp is None:
                self.QTapp = QApplication(sys.argv)
            else:
                print('QApplication allready exists: %s' % str(self.QTapp))
            if visTool == 3:
                slicerTool.slicerInterface(self.vtiFile)
            else:
                cutoffTool.cutoffInterface(self.vtiFile)
        else:
            raise NameError('No valid view option specified')

    def generateFullFlowfield(self):
        xLen = np.linspace(-250, max(self.output.rotLocX) + 1000, 200)
        yLen = np.linspace(-250, max(self.output.rotLocY) + 250, 200)
        zLen = np.linspace(0, 2*max(self.output.layout.locZ), 50)
        X, Y, Z = np.meshgrid(xLen, yLen, zLen, indexing='ij')
        U = florisCoreFunctions.windPlant.velAtLocations(X, Y, Z, self.output)
        return flowData(X, Y, Z, U)

    def generateHubFlowfield(self):
        D0 = self.output.layout.turbines[0].rotorDiameter
        xLen = np.linspace(-2*D0, max(self.output.rotLocX) + 15*D0, 200)
        yLen = np.linspace(-2*D0, max(self.output.rotLocY) + 2*D0, 200)
        zLen = self.output.layout.locZ[0]
        X, Y, Z = np.meshgrid(xLen, yLen, zLen, indexing='ij')
        U = florisCoreFunctions.windPlant.velAtLocations(X, Y, Z, self.output)
        return flowData(X, Y, Z, U)

    def generateSliceFlowfield(self, x, turbI):
        xLen = self.output.rotLocX[turbI] + x
        yLen = np.linspace(-250, max(self.output.rotLocY) + 250, 200)
        zLen = np.linspace(0, 2*max(self.output.layout.locZ), 50)
        X, Y, Z = np.meshgrid(xLen, yLen, zLen, indexing='ij')
        U = florisCoreFunctions.windPlant.velAtLocations(X, Y, Z, self.output)
        return flowData(X, Y, Z, U)

    def generateLidarFlowfield(self, turbI):
        # TODO: This might also need to be implemented
        #vis['LidarParams'] = 'LiDAR/Liss2Grid7x7widescreen_1d000_0d450_2d800_10000_1d600_NRELGE_Continuous_Parameter_JENExport.mat'
        #vis['LidarTrajectory'] = 'LiDAR/Liss2Grid7x7widescreen_1d000_0d450_2d800_10000_1d600_NRELGE_Continuous_JENExport.mat'

        XList = read_csv(self.output.prPath + '/LiDAR/xLoc.csv').values
        YList = read_csv(self.output.prPath + '/LiDAR/yLoc.csv').values
        ZList = read_csv(self.output.prPath + '/LiDAR/zLoc.csv').values

        # TODO, this does not represent xList yList and zList, fix it
        xLen = self.output.rotLocX[turbI] + np.unique(XList)
        yLen = self.output.rotLocY[turbI] + np.unique(YList)
        zLen = self.output.layout.locZ[turbI] + np.unique(ZList)
        X, Y, Z = np.meshgrid(xLen, yLen, zLen, indexing='ij')
        U = florisCoreFunctions.windPlant.velAtLocations(X, Y, Z, self.output)
        return flowData(X, Y, Z, U)
