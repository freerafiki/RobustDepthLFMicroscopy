function compileCPPFiles()

cd DEPTH/CPP/
mex tri_stereo_v2.cpp
mex bin_stereo_v3.cpp
mex compute_diff_image_aw.cpp
cd ../..

cd DENOISING/CPP/
mex jointBF_cpp.cpp
cd ../..
%cd SUPERPIXELS/SLIC_mex/
%mex slicmex.c
%mex slicomex.c
%cd ../..