function wmf_img = wmf(img, ref_img, ws)

    %% COMPUTE WEIGHTED MEDIAN FILTER (with non-integer weights) 
    tic
    cd DENOISING/CPP/
    mex wmf_cpp.cpp
    cd ../..
    wmf_img = wmf_cpp(double(img), double(ref_img), ws);
    toc
end