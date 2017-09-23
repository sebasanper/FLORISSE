# Notes on the FLORIS toolbox
-------------------------------

This FLORIS model is designed to provide a computationally inexpensive, controls-oriented model of the steady-state wake characteristics in a wind farm.  This can be used for online optimization.  This version of FLORIS implements a 3D version of the Jensen, original FLORIS (Gebraad et. al. 2016), and Gaussian wake model.  Literature on the Gaussian model can be found in the following papers:

1. Niayifar, A. and Port�-Agel, F.: A new 15 analytical model for wind farm power prediction, in: Journal of Physics: Conference Series, vol. 625, p. 012039, IOP Publishing, 2015.

2. Dilip, D. and Port�-Agel, F.: Wind Turbine Wake Mitigation through Blade Pitch Offset, Energies, 10, 757, 2017.

3. Abkar, M. and Port�-Agel, F.: Influence of atmospheric stability on wind-turbine wakes: A large-eddy simulation study, Physics of Fluids, 27, 035 104, 2015.

4. Bastankhah, M. and Port�-Agel, F.: A new analytical model for wind-turbine wakes, Renewable Energy, 70, 116�123, 2014.

5. Bastankhah, M. and Port�-Agel, 5 F.: Experimental and theoretical study of wind turbine wakes in yawed conditions, Journal of FluidMechanics, 806, 506�541, 2016.


System Requirements:
--------------------
Python 3

Numpy

Scipy 

Pandas

Matplotlib

Jupyter Notebooks

System Requirements for optional VTK visualization tools:
--------------------
vtk (Recomended to use Anaconda 64bit, install VTK using
     $ conda install -c clinicalgraphics vtk
     )

To Run:
--------
Open RunGAUSS_3D.ipynb - this notebook contains 5 examples of how to use this FLORIS code:

	1. Extract Velocity and Power information

	2. Visualize the wind farm in the horizontal plane using different wake models

	3. Wake steering optimization (objective of power maximization)

	4. Thrust optimization (objective of power maximization)

	5. Visualize Lidar module based on the University of Stuttgart scanning Lidar

	6. Extract velocities at specified points in the field


**Note: a detailed inputData dictionary is necessary to run the code.  This is detailed in RunGAUSS_3D.ipynb

 

