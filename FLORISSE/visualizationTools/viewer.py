# -*- coding: utf-8 -*-
import numpy as np
import sys
import os.path

import florisCoreFunctions.windPlant
import visualizationTools.MPLVisualizations as MPLViews
import visualizationTools.cutoffVisualizationTool as cutoffTool
import visualizationTools.slicerTool as slicerTool
from visualizationTools.writeVtiFlowfield import vtiFlowfieldWriter
from PyQt5.QtWidgets import QApplication


class viewer:
    def __init__(self, output):
        self.prepared = []
        self.vtiFile = []
        self.output = output

    def showView(self, visTool):
        if visTool not in self.prepared:
            print('prepare')

            if (visTool == 3 or visTool == 4) and not self.vtiFile:
                self.vtiFile = (self.output.rootPath +
                                '/visualizationData/currentFlowData.vti')
                vFW = vtiFlowfieldWriter(self.vtiFile)
                X, Y, Z, U = self.generateFullFlowfield()
                vFW.write(X, Y, Z, U)
                self.prepared.append([3, 4])

        if visTool == 0:
            MPLViews.visualizeHorizontal
        elif visTool == 1:
            MPLViews.visualizeCut
        elif visTool == 2:
            MPLViews.visualizeLidar
        elif visTool == 3 or visTool == 4:
            self.app = QApplication.instance()
            if self.app is None:
                self.app = QApplication(sys.argv)
            else:
                print('QApplication instance already exists: %s' % str(self.app))
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
        Z, Y, X = np.meshgrid(zLen, yLen, xLen, indexing='ij')
        U = florisCoreFunctions.windPlant.velAtLocations(X, Y, Z, self.output)
        return X, Y, Z, U
