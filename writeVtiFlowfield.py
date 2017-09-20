# -*- coding: utf-8 -*-
import numpy as np
import vtk

xVec = outputData['xLen']
yVec = outputData['yLen']
zVec =np.linspace(inputData['zLen'][0],inputData['zLen'][1],inputData['nSamplesZ'])
U = outputData['Ufield']


imageData = vtk.vtkImageData()
imageData.SetDimensions(inputData['nSamplesX'], inputData['nSamplesY'],
                        inputData['nSamplesZ'])
imageData.AllocateScalars(vtk.VTK_DOUBLE, 1)
imageData.SetSpacing(xVec[1]-xVec[0],yVec[1]-yVec[0],zVec[1]-zVec[0])
dims = imageData.GetDimensions()

# Fill every entry of the image data with "2.0"
for z in range(dims[2]):
    for y in range(dims[1]):
        for x in range(dims[0]):
            imageData.SetScalarComponentFromDouble(x, y, z, 0, U[z,y,x])
 
writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName("flowData.vti")
writer.SetInputData(imageData)
writer.Write()