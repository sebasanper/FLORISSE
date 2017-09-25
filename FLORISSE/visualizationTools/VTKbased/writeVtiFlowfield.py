# -*- coding: utf-8 -*-
import vtk


class vtiFlowfieldWriter:
    def __init__(self, filename):
        self.writer = vtk.vtkXMLImageDataWriter()
        self.writer.SetFileName(filename)

    def write(self, flowField):
        xVec = flowField.X[:, 0, 0]
        yVec = flowField.Y[0, :, 0]
        zVec = flowField.Z[0, 0, :]

        imageData = vtk.vtkImageData()
        imageData.SetDimensions(flowField.U.shape)
        imageData.AllocateScalars(vtk.VTK_DOUBLE, 1)
        imageData.SetSpacing(xVec[1]-xVec[0], yVec[1]-yVec[0], zVec[1]-zVec[0])
        dims = imageData.GetDimensions()

        # Fill every entry of the image data with "2.0"
        for x in range(dims[0]):
            for y in range(dims[1]):
                for z in range(dims[2]):
                    imageData.SetScalarComponentFromDouble(
                            x, y, z, 0, flowField.U[x, y, z])
        self.writer.SetInputData(imageData)
        self.writer.Write()
