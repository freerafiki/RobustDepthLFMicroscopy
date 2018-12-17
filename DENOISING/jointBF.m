function denorm_img = jointBF(img, ref_img, ws)

    %% COMPUTE WEIGHTED MEDIAN FILTER (with non-integer weights) 
    %tic
    %cd DENOISING/CPP/
    %mex jointBF_cpp.cpp
    %cd ../..
    grayref = rgb2gray(uint8(ref_img));
    maxval = max(max(img));
    img_norm = img .* (255/maxval); % bring to 255 for the weights
    filtered_img = jointBF_cpp(double(img_norm), double(grayref), ws);
    denorm_img = filtered_img ./ (255/maxval);
    %toc
end