%%%%%%%%%%%%%%%%%%%%%
%
% IDEA: compute three depth pairs (from the three directions)
%
% STEPS:
% 1) compute cost cubes in the three directions across the three scales
% 2) sum them up using frequency map - 6 cubes
% 3) compute the depth maps using MGM
% 4) apply left-right consistency check for eliminating outliers
% 5) fuse the three depths using failure prediciton map
%
% PROBLEMS:
% where do I use superpixels and how?
%
%%%%%%%%%%%%%%%%%%%%

function [cv, depthmap, conf] = hex_stereo_matching_v2(EIs, c_img, ws, alpha, d_min, d_max, matte, gamma1, gamma2, priors, use_failure_prediction)

    fprintf('\nComputing cost volume from correspondences!\n')
    % computes disparity from correspondences using multi view approach
    tic
       
    fprintf('First Epipolar Line..\n');
    %% 1) compute cost cubes
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
    imgR = rgb2gray(uint8(EIs(:,:,:,7)));
    [cvL2, cvR2] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 3);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL2_s1, cvR2_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 3);
    [cvL2_s2, cvR2_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 3);
    
    fprintf('Third Epipolar Line..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,2)));
    imgR = rgb2gray(uint8(EIs(:,:,:,6)));
    [cvL3, cvR3] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 2);
    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL3_s1, cvR3_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 2);
    [cvL3_s2, cvR3_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 2);
    
    %% 2) sum them up
    fprintf('Priors..\n');
    %% PREPARE PRIORS MAPS
    thresh = multithresh(priors,2);
    maps = imquantize(priors,thresh);
    high_freqs = maps == 3;
    med_freqs = maps == 2;
    low_freqs = maps == 1;

    % six cubes
    cvL1msc = zeros(size(cvL1,1), size(cvL1,2), size(cvL1,3));
    cvL2msc = zeros(size(cvL1,1), size(cvL1,2), size(cvL1,3));
    cvL3msc = zeros(size(cvL1,1), size(cvL1,2), size(cvL1,3));
    cvR1msc = zeros(size(cvL1,1), size(cvL1,2), size(cvL1,3));
    cvR2msc = zeros(size(cvL1,1), size(cvL1,2), size(cvL1,3));
    cvR3msc = zeros(size(cvL1,1), size(cvL1,2), size(cvL1,3));
    
    
    fprintf('Aggregating Multi-Scale Costs..\n');
    r = 5;
    eps = 0.001;
    guidedfilter_color_precompute(c_img, r, eps);
    for m = 1: size(cvL1, 3)

        % normal depth resolution scale 0
        slice_L1_s0 = cvL1(:,:,m);
        slice_L2_s0 = cvL2(:,:,m);
        slice_L3_s0 = cvL3(:,:,m);
        slice_R1_s0 = cvR1(:,:,m);
        slice_R2_s0 = cvR2(:,:,m);
        slice_R3_s0 = cvR3(:,:,m);

        % reduced scale 1
        m1 = min(max(1,int8(floor(m/2))),size(cvL1_s1,3));
        slice_L1_s1 = impyramid(cvL1_s1(:,:,m1),'expand');
        slice_L2_s1 = impyramid(cvL2_s1(:,:,m1),'expand');
        slice_L3_s1 = impyramid(cvL3_s1(:,:,m1),'expand');
        slice_R1_s1 = impyramid(cvR1_s1(:,:,m1),'expand');
        slice_R2_s1 = impyramid(cvR2_s1(:,:,m1),'expand');
        slice_R3_s1 = impyramid(cvR3_s1(:,:,m1),'expand');

        % reduced scale 2
        m2 = min(max(1, int8(round(m/4))), size(cvL1_s2,3));
        slice_L1_s2_tmp = impyramid(impyramid(cvL1_s2(:,:,m2),'expand'),'expand');
        slice_L2_s2_tmp = impyramid(impyramid(cvL2_s2(:,:,m2),'expand'),'expand');
        slice_L3_s2_tmp = impyramid(impyramid(cvL3_s2(:,:,m2),'expand'),'expand');
        slice_R1_s2_tmp = impyramid(impyramid(cvR1_s2(:,:,m2),'expand'),'expand');
        slice_R2_s2_tmp = impyramid(impyramid(cvR2_s2(:,:,m2),'expand'),'expand');
        slice_R3_s2_tmp = impyramid(impyramid(cvR3_s2(:,:,m2),'expand'),'expand');
        if size(slice_L1_s2_tmp,1) < size(slice_L1_s1,1)
            slice_L1_s2 = zeros(size(slice_L1_s1,2), size(slice_L1_s1,2));
            slice_L1_s2(2:size(slice_L1_s2,1)-1, 2:size(slice_L1_s2,2)-1) = slice_L1_s2_tmp;
            slice_L2_s2 = zeros(size(slice_L1_s1,2), size(slice_L1_s1,2));
            slice_L2_s2(2:size(slice_L2_s2,1)-1, 2:size(slice_L2_s2,2)-1) = slice_L2_s2_tmp;
            slice_L3_s2 = zeros(size(slice_L1_s1,2), size(slice_L1_s1,2));
            slice_L3_s2(2:size(slice_L3_s2,1)-1, 2:size(slice_L3_s2,2)-1) = slice_L3_s2_tmp;
            slice_R1_s2 = zeros(size(slice_L1_s1,2), size(slice_L1_s1,2));
            slice_R1_s2(2:size(slice_R1_s2,1)-1, 2:size(slice_R1_s2,2)-1) = slice_R1_s2_tmp;
            slice_R2_s2 = zeros(size(slice_L1_s1,2), size(slice_L1_s1,2));
            slice_R2_s2(2:size(slice_R2_s2,1)-1, 2:size(slice_R2_s2,2)-1) = slice_R2_s2_tmp;
            slice_R3_s2 = zeros(size(slice_L1_s1,2), size(slice_L1_s1,2));
            slice_R3_s2(2:size(slice_R3_s2,1)-1, 2:size(slice_R3_s2,2)-1) = slice_R3_s2_tmp;
        else
            slice_L1_s2 = slice_L1_s2_tmp;
            slice_L2_s2 = slice_L2_s2_tmp;
            slice_L3_s2 = slice_L3_s2_tmp;
            slice_R1_s2 = slice_R1_s2_tmp;
            slice_R2_s2 = slice_R2_s2_tmp;
            slice_R3_s2 = slice_R3_s2_tmp;
        end

        % sum them up
        cur_slice_L1 = high_freqs .* (slice_L1_s0 .* 0.6 + slice_L1_s1 .* 0.3 + slice_L1_s2 .* 0.1) + ...
            med_freqs .* (slice_L1_s0 .* 0.2 + slice_L1_s1 .* 0.6 + slice_L1_s2 .* 0.2) + ...
            low_freqs .* (slice_L1_s0 .* 0.1 + slice_L1_s1 .* 0.3 + slice_L1_s2 .* 0.6);
        
        cur_slice_L2 = high_freqs .* (slice_L2_s0 .* 0.6 + slice_L2_s1 .* 0.3 + slice_L2_s2 .* 0.1) + ...
            med_freqs .* (slice_L2_s0 .* 0.2 + slice_L2_s1 .* 0.6 + slice_L2_s2 .* 0.2) + ...
            low_freqs .* (slice_L2_s0 .* 0.1 + slice_L2_s1 .* 0.3 + slice_L2_s2 .* 0.6);
        
        cur_slice_L3 = high_freqs .* (slice_L3_s0 .* 0.6 + slice_L3_s1 .* 0.3 + slice_L3_s2 .* 0.1) + ...
            med_freqs .* (slice_L3_s0 .* 0.2 + slice_L3_s1 .* 0.6 + slice_L3_s2 .* 0.2) + ...
            low_freqs .* (slice_L3_s0 .* 0.1 + slice_L3_s1 .* 0.3 + slice_L3_s2 .* 0.6);
        
        cur_slice_R1 = high_freqs .* (slice_R1_s0 .* 0.6 + slice_R1_s1 .* 0.3 + slice_R1_s2 .* 0.1) + ...
            med_freqs .* (slice_R1_s0 .* 0.2 + slice_R1_s1 .* 0.6 + slice_R1_s2 .* 0.2) + ...
            low_freqs .* (slice_R1_s0 .* 0.1 + slice_R1_s1 .* 0.3 + slice_R1_s2 .* 0.6);
        
        cur_slice_R2 = high_freqs .* (slice_R2_s0 .* 0.6 + slice_R2_s1 .* 0.3 + slice_R2_s2 .* 0.1) + ...
            med_freqs .* (slice_R2_s0 .* 0.2 + slice_R2_s1 .* 0.6 + slice_R2_s2 .* 0.2) + ...
            low_freqs .* (slice_R2_s0 .* 0.1 + slice_R2_s1 .* 0.3 + slice_R2_s2 .* 0.6);
        
        cur_slice_R3 = high_freqs .* (slice_R3_s0 .* 0.6 + slice_R3_s1 .* 0.3 + slice_R3_s2 .* 0.1) + ...
            med_freqs .* (slice_R3_s0 .* 0.2 + slice_R3_s1 .* 0.6 + slice_L1_s2 .* 0.2) + ...
            low_freqs .* (slice_R3_s0 .* 0.1 + slice_R3_s1 .* 0.3 + slice_L1_s2 .* 0.6);
       	

        % apply matte and guided filtering
        cvL1msc(:,:,m) = min(max(0, guidedfilter_color_runfilter(cur_slice_L1 .* double(matte))), 1);
        cvL2msc(:,:,m) = min(max(0, guidedfilter_color_runfilter(cur_slice_L2 .* double(matte))), 1);
        cvL3msc(:,:,m) = min(max(0, guidedfilter_color_runfilter(cur_slice_L3 .* double(matte))), 1);
        cvR1msc(:,:,m) = min(max(0, guidedfilter_color_runfilter(cur_slice_R1 .* double(matte))), 1);
        cvR2msc(:,:,m) = min(max(0, guidedfilter_color_runfilter(cur_slice_R2 .* double(matte))), 1);
        cvR3msc(:,:,m) = min(max(0, guidedfilter_color_runfilter(cur_slice_R3 .* double(matte))), 1);

        
    end
    
    %% 3) compute depth map using MGM
    %NDIR = 8; % Number of directions
    %P1 = 0.1; % Penalty 1
    %P2 = 0.2; % Penalty 2
    %MGM = 4; % 1, 2 or 4 messages
    %VTYPE = 2; % TLP

    %w = create_edge_costs(c_img);

    %% MGMs
    %{
    addpath(genpath('/data1/palmieri/MGM'));
    cd /data1/palmieri/MGM
    depth_L1 = MGM_wrapper(cvL1msc, NDIR, P1, P2, MGM, VTYPE, w);
    depth_R1 = MGM_wrapper(cvR1msc, NDIR, P1, P2, MGM, VTYPE, w);
    depth_L2 = MGM_wrapper(cvL2msc, NDIR, P1, P2, MGM, VTYPE, w);
    depth_R2 = MGM_wrapper(cvR2msc, NDIR, P1, P2, MGM, VTYPE, w);
    depth_L3 = MGM_wrapper(cvL3msc, NDIR, P1, P2, MGM, VTYPE, w);
    depth_R3 = MGM_wrapper(cvR3msc, NDIR, P1, P2, MGM, VTYPE, w);
    cd /data1/palmieri/Valencia/Code/matlab_depth
    %}
    
    c_img = double(rgb2gray(uint8(c_img)));
    %% PREPARE FAILURE MAPS
    if use_failure_prediction 
        % create the 3 maps
        kernel_0_deg = [0.5 1 0 -1 -0.5; 1 2 0 -2 -1; 2 3 0 -3 -2; 1 2 0 -2 -1; 0.5 1 0 -1 -0.5];
        kernel_60_deg = [1 0 -0.5 -2 -1; 2 0 0 -2 -1; 2 1 0 -1 -2; 1 2 0 0 -2; 1 2 0.5 0 -1];
        kernel_120_deg = [1 2 0.5 0 -1; 1 2 0 0 -2; 2 1 0 -1 -2; 2 0 0 -2 -1; 1 0 -0.5 -2 -1];
        epi1_fail_map = conv2(c_img, kernel_0_deg', 'same');
        epi2_fail_map = conv2(c_img, kernel_60_deg, 'same');
        epi3_fail_map = conv2(c_img, kernel_120_deg, 'same');
        %c_img_s1 = impyramid(c_img, 'reduce');
        %c_img_s2 = impyramid(c_img_s1, 'reduce');
        %epi1_fail_map_s1 = conv2(c_img_s1, kernel_0_deg, 'same');
        %epi1_fail_map_s2 = conv2(c_img_s2, kernel_0_deg, 'same');
        %epi2_fail_map_s1 = conv2(c_img_s1, kernel_60_deg, 'same');
        %epi2_fail_map_s2 = conv2(c_img_s2, kernel_60_deg, 'same');
        %epi3_fail_map_s1 = conv2(c_img_s1, kernel_120_deg, 'same');
        %epi3_fail_map_s2 = conv2(c_img_s2, kernel_120_deg, 'same');
            
        max_epi1 = max(max(epi1_fail_map));
        max_epi2 = max(max(epi2_fail_map));
        max_epi3 = max(max(epi3_fail_map));
        max_val = max(max(max_epi1, max_epi2), max_epi3);
        quantiz_levels = 10;
        step = max_val / quantiz_levels;
        vec_1 = reshape(abs(epi1_fail_map), size(epi1_fail_map, 1) * size(epi1_fail_map, 2), 1);
        fail1 = Quantiz(vec_1, 0:step:max_val);
        fail1_map = reshape(fail1, [size(epi1_fail_map,1), size(epi1_fail_map,2)]) ;
        
        vec_2 = reshape(abs(epi2_fail_map), size(epi2_fail_map, 1) * size(epi2_fail_map, 2), 1);
        fail2 = Quantiz(vec_2, 0:step:max_val);
        fail2_map = reshape(fail2, [size(epi2_fail_map,1), size(epi2_fail_map,2)]) ;
        
        vec_3 = reshape(abs(epi3_fail_map), size(epi3_fail_map, 1) * size(epi3_fail_map, 2), 1);
        fail3 = Quantiz(vec_3, 0:step:max_val);
        fail3_map = reshape(fail3, [size(epi3_fail_map,1), size(epi3_fail_map,2)]) ;
    end
    
    %% what if I median the costs and then extract?
    cur_cubes = zeros(size(cvL1msc, 1), size(cvL1msc, 2), 6);
    new_cv = zeros(size(cvL1msc, 1), size(cvL1msc, 2), size(cvL1msc, 3));
    var_cv = zeros(size(cvL1msc, 1), size(cvL1msc, 2), size(cvL1msc, 3));
    med_cv = zeros(size(cvL1msc, 1), size(cvL1msc, 2), size(cvL1msc, 3));
    cv = zeros(size(cvL1msc, 1), size(cvL1msc, 2), size(cvL1msc, 3));
    
    for k = 1: size(cvL1msc,3)
        
        cur_cubes(:,:,1) = cvL1msc(:,:,k);
        cur_cubes(:,:,2) = cvL2msc(:,:,k);
        cur_cubes(:,:,3) = cvL3msc(:,:,k);
        cur_cubes(:,:,4) = cvR1msc(:,:,k);
        cur_cubes(:,:,5) = cvR2msc(:,:,k);
        cur_cubes(:,:,6) = cvR3msc(:,:,k);
        
        if use_failure_prediction
            total_fail = fail1_map + fail2_map + fail3_map;
            cv(:,:,k) = (cvL1msc(:,:,k) .* fail1_map./total_fail + cvR1msc(:,:,k) .* fail1_map./total_fail) ./ 2 + ...
                (cvL2msc(:,:,k) .* fail2_map./total_fail + cvR2msc(:,:,k) .* fail2_map./total_fail) ./ 2 + ...
                (cvL3msc(:,:,k) .* fail3_map./total_fail + cvR3msc(:,:,k) .* fail3_map./total_fail) ./ 2;
        else
            cv(:,:,k) = mean(cur_cubes, 3);
        end
        
        var_cv(:,:,k) = var(cur_cubes, 0, 3);
        med_cv(:,:,k) = median(cur_cubes, 3);
        %cv(:,:,k) = mean(cur_cubes, 3);
        new_cv(:,:,k) = imabsdiff(med_cv(:,:,k), cv(:,:,k));
    end

    %% WTAs
    [~, depth_L1_wta] = min(cvL1msc, [], 3);
    [~, depth_R1_wta] = min(cvR1msc, [], 3);
    [~, depth_L2_wta] = min(cvL2msc, [], 3);
    [~, depth_R2_wta] = min(cvR2msc, [], 3);
    [~, depth_L3_wta] = min(cvL3msc, [], 3);
    [~, depth_R3_wta] = min(cvR3msc, [], 3);
    
    %% apply matte?
    apply_matte = true;
    if apply_matte 
        
        depth_L1_wta = depth_L1_wta .* double(matte);
        depth_R1_wta = depth_R1_wta .* double(matte);
        depth_L2_wta = depth_L2_wta .* double(matte);
        depth_R2_wta = depth_R2_wta .* double(matte);
        depth_L3_wta = depth_L3_wta .* double(matte);
        depth_R3_wta = depth_R3_wta .* double(matte);
    
    end

    %% 4) Left-Right Consistency Check 
    %structure with depths
    depths = zeros(size(depth_L1_wta, 1), size(depth_L1_wta,2), 6);
    depths(:,:,1) = depth_L1_wta;
    depths(:,:,2) = depth_L2_wta;
    depths(:,:,3) = depth_L3_wta;
    depths(:,:,4) = depth_R1_wta;
    depths(:,:,5) = depth_R2_wta;
    depths(:,:,6) = depth_R3_wta;
    
    depthmap = median(depths, 3);
    test = mean(depths, 3);
    conff = imabsdiff(depthmap, test);
    conf = (max(max(conff)) - conff) .* double(matte);
    
    %% normalize the confidence
    conf = conf ./ (max(max(conf)) + 1);
    
    %% invert the depth
    depthmap = (size(cvL1,3) - depthmap) .* double(matte);

        
    %{
    [~, vardepth] = min(var_cv, [], 3);
    [~, expdepth] = min(new_cv, [], 3);
    [~, meddepth] = min(med_cv, [], 3);
    [~, meandepth] = min(cv, [], 3);
    
    subplot(221); imagesc(depthmap);
    subplot(222); imagesc(expdepth);
    subplot(223); imagesc(meddepth);
    subplot(224); imagesc(meandepth);
    
   
    %% showing off
    subplot(3,2,1)
    imagesc(depth_L1)
    subplot(3,2,2)
    imagesc(depth_R1)
    subplot(3,2,3)
    imagesc(depth_L2)
    subplot(3,2,4)
    imagesc(depth_R2)
    subplot(3,2,5)
    imagesc(depth_L3)
    subplot(3,2,6)
    imagesc(depth_R3)
    
    figure()
    subplot(3,2,1)
    imagesc(depth_L1_wta)
    subplot(3,2,2)
    imagesc(depth_R1_wta)
    subplot(3,2,3)
    imagesc(depth_L2_wta)
    subplot(3,2,4)
    imagesc(depth_R2_wta)
    subplot(3,2,5)
    imagesc(depth_L3_wta)
    subplot(3,2,6)
    imagesc(depth_R3_wta)
    
    figure()
    imagesc(final_depth)
    %}
    fprintf('Done!                       ');
    toc
    
    