# Robust Depth Estimation for Light Field Microscopy

This is a Matlab implementation of the procedure described in the paper:

*Robust Depth Estimation for Light Field Microscopy*, L. Palmieri, G. Scrofani, N. Incardona, G. Saavedra, M. Martínez-Corral and R. Koch, Sensors 2019

Paper has been accepted and will be soon uploaded for more information.

It combines vision cues (correspondence, defocus) to estimate accurate and robust depth maps from light field microscopy images. Images are acquired with FiMic microscope or created with Blender.

### Microscope Imagery
The FiMic microscope allows to capture light field in a single shot.
For further information regarding FiMic Fourier Integral Microscope, refer to the paper [FIMic: design for ultimate 3D-integral microscopy of in-vivo biological samples](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5772586/)

### Implementation
About implementation, the code is mainly in Matlab and it makes use of some MEX files (C and C++) to speed up some computations. It uses other toolboxes, but they are already present in the github folders (see References for more details).

#### Dependencies
The dependencies are a C++ compiler for the MEX code (if you are using gcc, supported version is 4.7; however until 6.x is working and Matlab just gives a warning, from 7.0 it may cause some problems. If on Linux systems, use of alternatives can help about this [read more here](https://askubuntu.com/questions/26498/how-to-choose-the-default-gcc-and-g-version))

#### References
It takes inspiration and parts of code from other approaches, who have been cited in the paper and we would like to thank for sharing their code.

[1] JEON2015: [Open Source Code](https://drive.google.com/file/d/0B2553ggh3QTcS01zU0RjOG5FTjQ/view)
 -- [Accurate Depth Map Estimation from a Lenslet Light Field Camera, Jeon et al., 2015](https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Jeon_Accurate_Depth_Map_2015_CVPR_paper.pdf)

[2] MGM2015: [Project page](http://www.bmva.org/bmvc/2015/papers/paper090/index.html) -- [Open Source Code](http://www.bmva.org/bmvc/2015/papers/paper090/index.html) -- 
[MGM: A Significantly More Global Matching for Stereovision, Facciolo et al., 2015](http://www.bmva.org/bmvc/2015/papers/paper090/paper090.pdf)

[3] SFF2016: [Open Source Code](https://sites.google.com/view/cvia/downloads) -- [Analysis of focus measure operators for shape-from-focus, Pertuz et al., 2013](https://www.sciencedirect.com/science/article/pii/S0031320312004736?via%3Dihub)

[4] TAO2013: [Project Page with Open Source Code](http://graphics.berkeley.edu/papers/Tao-DFC-2013-12/) -- [Depth from Combining Defocus and Correspondence Using Light-Field Cameras, Tao et al., 2013](http://graphics.berkeley.edu/papers/Tao-DFC-2013-12/Tao-DFC-2013-12.pdf)


For further information, please contact me at: lpa@informatik.uni-kiel.de
