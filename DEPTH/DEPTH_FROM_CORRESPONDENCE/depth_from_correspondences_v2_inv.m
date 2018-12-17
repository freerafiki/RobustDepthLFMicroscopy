function [cvMS, cvMSGF] = depth_from_correspondences_v2(EIs, c_img, ws, alpha, d_min, d_max, matte, gamma_s1, gamma_s2)

    fprintf('\nComputing depth from correspondences!\n')
    % computes disparity from correspondences using multi view approach
    tic
    num_of_imgs = size(EIs,4);
    big_matrix = zeros(size(c_img,1), size(c_img,2), num_of_imgs);
    hws = floor(ws/2);
    selpix = [150, 350];
    selpix_cost = zeros(num_of_imgs,1);
    selpix_real_cost = zeros(num_of_imgs,1);
    for i=1:num_of_imgs
        gray_img = rgb2gray(uint8(EIs(:,:,:,i)));
        big_matrix(:,:,i) = gray_img;
    end
    gray_c_img = rgb2gray(uint8(c_img));
    disparities = linspace(d_min, d_max, d_max - d_min + 1);
    
    %imagesc(uint8(imgs{1,i}))
    %figure()
    %imagesc(rgb2gray(uint8(c_img)))
    cd DEPTH/CPP/
    %mex multi_view_CAD_v2.cpp
    mex tri_stereo.cpp
    cd ../..
    %cost_volume = multi_view_CAD_v2(double(big_matrix), (d_min), (d_max), (ws), (alpha));
    %[cost_volume_SAD, cost_volume_CENSUS] = multi_view_CAD(double(gray_c_img), double(big_matrix), double(disparities), (d_min), (d_max), (ws));
        
    fprintf('First Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,3)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,5)));
    [cvL1, cvR1] = tri_stereo(double(imgR), double(imgC), double(imgL), (d_min), (d_max), (ws), (alpha), 1);

    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL1_s1, cvR1_s1] = tri_stereo(double(imgR_s1), double(imgC_s1), double(imgL_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 1);
    [cvL1_s2, cvR1_s2] = tri_stereo(double(imgR_s2), double(imgC_s2), double(imgL_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 1);

    fprintf('Second Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,1)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,7)));
    [cvL2, cvR2] = tri_stereo(double(imgR), double(imgC), double(imgL), (d_min), (d_max), (ws), (alpha), 2);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL2_s1, cvR2_s1] = tri_stereo(double(imgR_s1), double(imgC_s1), double(imgL_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 2);
    [cvL2_s2, cvR2_s2] = tri_stereo(double(imgR_s2), double(imgC_s2), double(imgL_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 2);
    
    fprintf('Third Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,2)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,6)));
    [cvL3, cvR3] = tri_stereo(double(imgR), double(imgC), double(imgL), (d_min), (d_max), (ws), (alpha), 3);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL3_s1, cvR3_s1] = tri_stereo(double(imgR_s1), double(imgC_s1), double(imgL_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 3);
    [cvL3_s2, cvR3_s2] = tri_stereo(double(imgR_s2), double(imgC_s2), double(imgL_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 3);
    
    fprintf('Summing up contributions using multi-scale..\n');

    cvMS = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
    cvMSGF = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
    %gamma_s1 = 0.6;
    %gamma_s2 = 0.3;
    penalty = 0.1;    
    THRESH_DIFF = 0.01;
    for m = 1: size(cvL1, 3)
        % normal depth resolution scale 0
        LR1 = (cvL1(:,:,m) + cvR1(:,:,m)) / 2;
        LR2 = (cvL2(:,:,m) + cvR2(:,:,m)) / 2;
        LR3 = (cvL3(:,:,m) + cvR3(:,:,m)) / 2;
        different_pixels = ( (abs(LR1-LR2)>THRESH_DIFF) .* (abs(LR1-LR3)>THRESH_DIFF) .* (abs(LR1-LR3)>THRESH_DIFF) );
        slice_s0 = ( LR1 + LR2 + LR3 + penalty .* different_pixels ) / 3.0;
        % reduced scale 1
        m1 = min(max(1,int8(floor(m/2))),size(cvL1_s1,3));
        LR1_s1 = (cvL1_s1(:,:,m1) + cvR1_s1(:,:,m1)) / 2;
        LR2_s1 = (cvL2_s1(:,:,m1) + cvR2_s1(:,:,m1)) / 2;
        LR3_s1 = (cvL3_s1(:,:,m1) + cvR3_s1(:,:,m1)) / 2;
        different_pixels_s1 = ( (abs(LR1_s1-LR2_s1)>THRESH_DIFF) .* (abs(LR1_s1-LR3_s1)>THRESH_DIFF) .* (abs(LR1_s1-LR3_s1)>THRESH_DIFF) );
        slice_s1_LR = (LR1_s1 + LR2_s1 + LR3_s1 + penalty .* different_pixels_s1) / 3.0;
        slice_s1 = impyramid(slice_s1_LR,'expand');
        % reduced scale 2
        m2 = min(max(1, int8(round(m/4))), size(cvL1_s2,3));
        LR1_s2 = (cvL1_s2(:,:,m2) + cvR1_s2(:,:,m2)) / 2;
        LR2_s2 = (cvL2_s2(:,:,m2) + cvR2_s2(:,:,m2)) / 2;
        LR3_s2 = (cvL3_s2(:,:,m2) + cvR3_s2(:,:,m2)) / 2;
        different_pixels_s2 = ( (abs(LR1_s2-LR2_s2)>THRESH_DIFF) .* (abs(LR1_s2-LR3_s2)>THRESH_DIFF) .* (abs(LR1_s2-LR3_s2)>THRESH_DIFF) );
        slice_s2_LR = (LR1_s2 + LR2_s2 + LR3_s2 + penalty .* different_pixels_s2) / 3.0;
        slice_s2_tmp = impyramid(impyramid(slice_s2_LR,'expand'),'expand');
        if size(slice_s2_tmp,1) < size(slice_s1,1)
            slice_s2 = zeros(size(slice_s1,2), size(slice_s1,2));
            slice_s2(2:size(slice_s2,1)-1, 2:size(slice_s2,2)-1) = slice_s2_tmp;
        else
            slice_s2 = slice_s2_tmp;
        end
        % sum them up
        cvMS(:,:,m) = slice_s0 + gamma_s1.*slice_s1 + gamma_s2.*slice_s2;
        % filter with GF
        cvMSGF(:,:,m) = imguidedfilter(cvMS(:,:,m));
    end
    
    
    fprintf('Done!                       ');
    toc
end
