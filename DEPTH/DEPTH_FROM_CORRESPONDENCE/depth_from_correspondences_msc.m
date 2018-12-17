function cvMS = depth_from_correspondences_msc(EIs, c_img, ws, alpha, d_min, d_max, matte, gamma1, gamma2, priors, use_failure_prediction)

    fprintf('\nComputing cost volume from correspondences!\n')
    % computes disparity from correspondences using multi view approach
    tic
       
    fprintf('First Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,5)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,3)));
    [cvL1, cvR1] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 1);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL1_s1, cvR1_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 1);
    [cvL1_s2, cvR1_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 1);

    fprintf('Second Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,1)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,7)));
    [cvL2, cvR2] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 3);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL2_s1, cvR2_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 3);
    [cvL2_s2, cvR2_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 3);
    
    fprintf('Third Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,2)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,6)));
    [cvL3, cvR3] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 2);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL3_s1, cvR3_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 2);
    [cvL3_s2, cvR3_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 2);
    
    fprintf('Summing up contributions using multi-scale..\n');

    cvMS = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
    c_img = double(rgb2gray(uint8(c_img)));
    %% PREPARE FAILURE MAPS
    if use_failure_prediction 
        % create the 3 maps
        kernel_0_deg = [0.5 1 0 -1 -0.5; 1 2 0 -2 -1; 2 3 0 -3 -2; 1 2 0 -2 -1; 0.5 1 0 -1 -0.5];
        kernel_60_deg = [1 0 0.5 2 1; 2 0 0 2 2; 2 1 0 1 2; 2 2 0 0 2; 1 2 0.5 0 1];
        kernel_120_deg = [1 2 0.5 0 1; 2 2 0 0 2; 2 1 0 1 2; 2 0 0 2 2; 1 0 0.5 2 1];
        epi1_fail_map = conv2(c_img, kernel_0_deg, 'same');
        epi2_fail_map = conv2(c_img, kernel_60_deg, 'same');
        epi3_fail_map = conv2(c_img, kernel_120_deg, 'same');
        c_img_s1 = impyramid(c_img, 'reduce');
        c_img_s2 = impyramid(c_img_s1, 'reduce');
        epi1_fail_map_s1 = conv2(c_img_s1, kernel_0_deg, 'same');
        epi1_fail_map_s2 = conv2(c_img_s2, kernel_0_deg, 'same');
        epi2_fail_map_s1 = conv2(c_img_s1, kernel_60_deg, 'same');
        epi2_fail_map_s2 = conv2(c_img_s2, kernel_60_deg, 'same');
        epi3_fail_map_s1 = conv2(c_img_s1, kernel_120_deg, 'same');
        epi3_fail_map_s2 = conv2(c_img_s2, kernel_120_deg, 'same');
    end
    
    %% PREPARE PRIORS MAPS
    thresh = multithresh(priors,2);
    maps = imquantize(priors,thresh);
    high_freqs = maps == 3;
    med_freqs = maps == 2;
    low_freqs = maps == 1;
    
    penalty = 0.1;    
    THRESH_DIFF = 0.01;
    for m = 1: size(cvL1, 3)
        % normal depth resolution scale 0
        LR1 = (cvL1(:,:,m) + cvR1(:,:,m)) / 2;
        LR2 = (cvL2(:,:,m) + cvR2(:,:,m)) / 2;
        LR3 = (cvL3(:,:,m) + cvR3(:,:,m)) / 2;
        if use_failure_prediction
            slice_s0 = ( LR1 .* epi1_fail_map + LR2 .* epi2_fail_map + LR3 .* epi3_fail_map ) / 3.0;
        else
            different_pixels = ( (abs(LR1-LR2)>THRESH_DIFF) .* (abs(LR1-LR3)>THRESH_DIFF) .* (abs(LR1-LR3)>THRESH_DIFF) );
            slice_s0 = ( LR1 + LR2 + LR3 + penalty .* different_pixels ) / 3.0;
        end
        % reduced scale 1
        m1 = min(max(1,int8(floor(m/2))),size(cvL1_s1,3));
        LR1_s1 = (cvL1_s1(:,:,m1) + cvR1_s1(:,:,m1)) / 2;
        LR2_s1 = (cvL2_s1(:,:,m1) + cvR2_s1(:,:,m1)) / 2;
        LR3_s1 = (cvL3_s1(:,:,m1) + cvR3_s1(:,:,m1)) / 2;
        if use_failure_prediction
            slice_s1_LR = ( LR1_s1 .* epi1_fail_map_s1 + LR2_s1 .* epi2_fail_map_s1 + LR3_s1 .* epi3_fail_map_s1 ) / 3.0;
        else
            different_pixels_s1 = ( (abs(LR1_s1-LR2_s1)>THRESH_DIFF) .* (abs(LR1_s1-LR3_s1)>THRESH_DIFF) .* (abs(LR1_s1-LR3_s1)>THRESH_DIFF) );
            slice_s1_LR = (LR1_s1 + LR2_s1 + LR3_s1 + penalty .* different_pixels_s1) / 3.0;
        end
        slice_s1 = impyramid(slice_s1_LR,'expand');
        % reduced scale 2
        m2 = min(max(1, int8(round(m/4))), size(cvL1_s2,3));
        LR1_s2 = (cvL1_s2(:,:,m2) + cvR1_s2(:,:,m2)) / 2;
        LR2_s2 = (cvL2_s2(:,:,m2) + cvR2_s2(:,:,m2)) / 2;
        LR3_s2 = (cvL3_s2(:,:,m2) + cvR3_s2(:,:,m2)) / 2;
        if use_failure_prediction
            slice_s2_LR = ( LR1_s2 .* epi1_fail_map_s2 + LR2_s2 .* epi2_fail_map_s2 + LR3_s2 .* epi3_fail_map_s2 ) / 3.0;
        else
            different_pixels_s2 = ( (abs(LR1_s2-LR2_s2)>THRESH_DIFF) .* (abs(LR1_s2-LR3_s2)>THRESH_DIFF) .* (abs(LR1_s2-LR3_s2)>THRESH_DIFF) );
            slice_s2_LR = (LR1_s2 + LR2_s2 + LR3_s2 + penalty .* different_pixels_s2) / 3.0;
        end
        slice_s2_tmp = impyramid(impyramid(slice_s2_LR,'expand'),'expand');
        if size(slice_s2_tmp,1) < size(slice_s1,1)
            slice_s2 = zeros(size(slice_s1,2), size(slice_s1,2));
            slice_s2(2:size(slice_s2,1)-1, 2:size(slice_s2,2)-1) = slice_s2_tmp;
        else
            slice_s2 = slice_s2_tmp;
        end
        % sum them up
       	cur_slice = slice_s0 .* (gamma1 .* high_freqs) + slice_s0 .* (gamma2 .* med_freqs) + ...
            slice_s1 .* (gamma2 .* (high_freqs + low_freqs)) + slice_s1 .* (gamma1 .* med_freqs) + ...
            slice_s2 .* (gamma2 .* med_freqs) + slice_s0 .* (gamma1 .* low_freqs);
        cur_slice = cur_slice .* double(matte);
        % filter with GF
        cvMS(:,:,m) = imguidedfilter(cur_slice);

    end
    
    
    fprintf('Done!                       ');
    toc
end
