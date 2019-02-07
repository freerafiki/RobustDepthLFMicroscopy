# Robust Depth Estimation for Light Field Microscopy

This is a Matlab implementation of the procedure described in the paper:

[*Robust Depth Estimation for Light Field Microscopy*, L. Palmieri, G. Scrofani, N. Incardona, G. Saavedra, M. Martínez-Corral and R. Koch, Sensors 2019, 19, 500.](https://www.mdpi.com/1424-8220/19/3/500)



### Algorithm

![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/BlockDiagramv2.png)

Starting from a FiMic image, it creates a focal stack and separates the perspective views. Then two cost volumes are built, and refined using multi-scale and superpixels contributions.
Finally, using energy minimization the two volumes are combined and the final depth map is extracted.

#### Some Results
Cotton fibers           | Cotton fibers         | Zebrafish             | Chip                  | Solderings            |
|:---------------------:|:---------------------:|:---------------------:|:---------------------:|:---------------------:|
| ![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/1_-_img.png)| ![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/2_-_img.png)| ![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/zebrafish.jpg)| ![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/chip.png)| ![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/EI_4.png)|
|![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/1_-_depth.png)| ![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/2_-_depth.png)|![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/prop4simatting.png)|![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/imgGF_col_depth.png)|![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/soldingLUCA.png)|


### Microscope Imagery
The FiMic microscope allows to capture light field in a single shot.

![](https://github.com/PlenopticToolbox/RobustDepthLFMicroscopy/blob/master/readmeimages/fimic.jpg)

For further information regarding FiMic Fourier Integral Microscope, refer to the paper [FIMic: design for ultimate 3D-integral microscopy of in-vivo biological samples](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5772586/)

### Implementation
About implementation, the code is mainly in Matlab and it makes use of some MEX files (C and C++) to speed up some computations. It uses other toolboxes, but they are already present in the github folders (see References for more details).

#### Dependencies
The dependencies are a C++ compiler for the MEX code (if you are using gcc, supported version is 4.7; however until 6.x is working and Matlab just gives a warning, from 7.0 it may cause some problems. If on Linux systems, use of alternatives can help about this [read more here](https://askubuntu.com/questions/26498/how-to-choose-the-default-gcc-and-g-version))

#### Main References
It takes inspiration and parts of code from other approaches, who have been cited in the paper and we would like to thank for sharing their code.

[1] JEON2015: [Open Source Code](https://drive.google.com/file/d/0B2553ggh3QTcS01zU0RjOG5FTjQ/view)
 -- [Accurate Depth Map Estimation from a Lenslet Light Field Camera, Jeon et al., 2015](https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Jeon_Accurate_Depth_Map_2015_CVPR_paper.pdf)

[2] MGM2015: [Project page](http://www.bmva.org/bmvc/2015/papers/paper090/index.html) -- [Open Source Code](http://www.bmva.org/bmvc/2015/papers/paper090/index.html) -- 
[MGM: A Significantly More Global Matching for Stereovision, Facciolo et al., 2015](http://www.bmva.org/bmvc/2015/papers/paper090/paper090.pdf)

[3] SFF2016: [Open Source Code](https://sites.google.com/view/cvia/downloads) -- [Analysis of focus measure operators for shape-from-focus, Pertuz et al., 2013](https://www.sciencedirect.com/science/article/pii/S0031320312004736?via%3Dihub)

[4] TAO2013: [Project Page with Open Source Code](http://graphics.berkeley.edu/papers/Tao-DFC-2013-12/) -- [Depth from Combining Defocus and Correspondence Using Light-Field Cameras, Tao et al., 2013](http://graphics.berkeley.edu/papers/Tao-DFC-2013-12/Tao-DFC-2013-12.pdf)


For further information, please contact me at: lpa@informatik.uni-kiel.de
