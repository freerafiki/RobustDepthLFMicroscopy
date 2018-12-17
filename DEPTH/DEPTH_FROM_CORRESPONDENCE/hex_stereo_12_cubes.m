%%%%%%%%%%%%%%%%%%%%%
%
% IDEA: compute three depth pairs (from the three directions)
%
% STEPS:
% 1) compute cost cubes in the three directions across the three scales
% 2) sum them up using frequency map - 12 cubes
% 3) compute the depth maps using MGM
% 4) apply left-right consistency check for eliminating outliers
% 5) fuse the three depths using failure prediciton map
%
% PROBLEMS:
% where do I use superpixels and how?
%
%%%%%%%%%%%%%%%%%%%%

function cvMS = depth_from_correspondences_msc(EIs, c_img, ws, alpha, d_min, d_max, matte, gamma1, gamma2, priors, use_failure_prediction)

    fprintf('\nComputing cost volume from correspondences!\n')
    % computes disparity from correspondences using multi view approach
    tic
       
    fprintf('First Epipolar Line..\n');
    %% 1) compute cost cubes
    imgL = rgb2gray(uint8(EIs(:,:,:,5)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,3)));
    [LC1, CL1, CR1, RC1] = tri_stereo_LR(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 1);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [LC1_s1, CL1_s1, CR1_s1, RC1_s1] = tri_stereo_LR(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 1);
    [LC1_s2, CL1_s2, CR1_s2, RC1_s2] = tri_stereo_LR(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 1);

    fprintf('Second Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,7)));
    imgR = rgb2gray(uint8(EIs(:,:,:,1)));
    [LC2, CL2, CR2, RC2] = tri_stereo_LR(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 2);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [LC1_s1, CL1_s1, CR1_s1, RC1_s1] = tri_stereo_LR(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 2);
    [LC1_s2, CL1_s2, CR1_s2, RC1_s2] = tri_stereo_LR(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 2);
    
    fprintf('Third Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,6)));
    imgR = rgb2gray(uint8(EIs(:,:,:,2)));
    [LC3, CL3, CR3, RC3] = tri_stereo_LR(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 3);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [LC1_s1, CL1_s1, CR1_s1, RC1_s1] = tri_stereo_LR(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 3);
    [LC1_s2, CL1_s2, CR1_s2, RC1_s2] = tri_stereo_LR(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 3);
    
    %% 2) sum them up using frequency map
    %% PREPARE PRIORS MAPS
    thresh = multithresh(priors,2);
    maps = imquantize(priors,thresh);
    high_freqs = maps == 3;
    med_freqs = maps == 2;
    low_freqs = maps == 1;
    
    %% WE NEED 12 CUBES??
    
    penalty = 0.1;    
    THRESH_DIFF = 0.01;
    for m = 1: size(cvL1, 3)
        
        % reduced scale 1
        m1 = min(max(1,int8(floor(m/2))),size(cvL1_s1,3));
        
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