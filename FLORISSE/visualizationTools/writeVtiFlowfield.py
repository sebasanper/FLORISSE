# -*- coding: utf-8 -*-
import vtk


class vtiFlowfieldWriter:
    def __init__(self, filename):
        self.writer = vtk.vtkXMLImageDataWriter()
        self.writer.SetFileName(filename)

    def write(self, X, Y, Z, U):
        xVec = X[0, 0, :]
        yVec = Y[0, :, 0]
        zVec = Z[:, 0, 0]

        imageData = vtk.vtkImageData()
        imageData.SetDimensions(U.shape[2], U.shape[1], U.shape[0])
        imageData.AllocateScalars(vtk.VTK_DOUBLE, 1)
        imageData.SetSpacing(xVec[1]-xVec[0], yVec[1]-yVec[0], zVec[1]-zVec[0])
        dims = imageData.GetDimensions()

        # Fill every entry of the image data with "2.0"
        for z in range(dims[2]):
            for y in range(dims[1]):
                for x in range(dims[0]):
                    imageData.SetScalarComponentFromDouble(x, y, z, 0,
                                                           U[z, y, x])
        self.writer.SetInputData(imageData)
        self.writer.Write()
