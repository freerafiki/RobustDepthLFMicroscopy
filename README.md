# Robust Depth Estimation for Light Field Microscopy

This is a Matlab implementation of the procedure described in the paper:

*Robust Depth Estimation for Light Field Microscopy*, L. Palmieri, G. Scrofani, N. Incardona, G. Saavedra, M. Martínez-Corral and R. Koch, Sensors 2019

Paper has been accepted and will be soon uploaded for more information.

It combines vision cues (correspondence, defocus) to estimate accurate and robust depth maps from light field microscopy images. Images are acquired with FiMic microscope or created with Blender.

### Microscope Imagery
The FiMic microscope allows to capture light field in a single shot.
For further information regarding FiMic Fourier Integral Microscope, refer to the paper [FIMic: design for ultimate 3D-integral microscopy of in-vivo biological samples](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5772586/)

### Implementation
About implementation, the code is mainly in Matlab and it makes use of some MEX files (C and C++) to speed up some computations. It uses other toolboxes, but they are already present in the github folders.

#### Dependencies
The dependencies are a C++ compiler for the MEX code (if you are using gcc, supported version is 4.7; however until 6.x is working and Matlab just gives a warning, from 7.0 it may cause some problems. If on Linux systems, use of alternatives can help about this [read more here](https://askubuntu.com/questions/26498/how-to-choose-the-default-gcc-and-g-version))

For further information, please contact me at: lpa@informatik.uni-kiel.de
