function cost_volume = depth_from_correspondences(EIs, c_img, ws, alpha, d_min, d_max, matte)

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
    imgL = rgb2gray(uint8(EIs(:,:,:,3)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,5)));
    [cvL1, cvR1] = tri_stereo(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 1);

    imgL_s1 = impyramid(imgL, 'reduce');
    imgL_s2 = impyramid(imgL_s1, 'reduce');
    imgC_s1 = impyramid(imgC, 'reduce');
    imgC_s2 = impyramid(imgC_s1, 'reduce');
    imgR_s1 = impyramid(imgR, 'reduce');
    imgR_s2 = impyramid(imgR_s1, 'reduce');
    [cvL1_s1, cvR1_s1] = tri_stereo(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 1);
    [cvL1_s2, cvR1_s2] = tri_stereo(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 1);
    
    cvMS = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
    cvMSGF = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
    gamma_1 = 0.6;
    gamma_2 = 0.3;
    theta = 1;    
    for m = 1: size(cvL1, 3)
        % normal depth resolution scale 0
        slice_s0 = (cvL1(:,:,m) + cvR1(:,:,m)) / 2 + theta.*(abs(cvL1(:,:,m)-cvR1(:,:,m)));
        m1 = int8(round(m/2));
        slice_s1_LR = (cvL1_s1(:,:,m1) + cvR1_s1(:,:,m1)) / 2 + theta.*(abs(cvL1_s1(:,:,m1)-cvR1_s1(:,:,m1)));
        slice_s1 = impyramid(slice_s1_LR,'expand');
        m2 = max(1, int8(round(m/4)));
        slice_s2_LR = (cvL1_s2(:,:,m2) + cvR1_s2(:,:,m2)) / 2 + theta.*(abs(cvL1_s2(:,:,m2)-cvR1_s2(:,:,m2)));
        slice_s2_tmp = impyramid(impyramid(slice_s2_LR,'expand'),'expand');
        slice_s2 = zeros(size(slice_s1,2), size(slice_s1,2));
        slice_s2(2:size(slice_s2,1)-1, 2:size(slice_s2,2)-1) = slice_s2_tmp;
        cvMS(:,:,m) = slice_s0 + gamma_1.*slice_s1 + gamma_2.*slice_s2;
        cvMSGF(:,:,m) = imguidedfilter(cvMS(:,:,m));
    end
    
    [d1, dLR1] = min((cvL1+cvR1)/2, [], 3);
    [d1, dL1_s1] = min(cvL1_s1, [], 3);
    [d1, dL1_s2] = min(cvL1_s2, [], 3);
    [dm, dMS] = min(cvMS, [], 3);
    [dg, dMSGF] = min(cvMSGF, [], 3);
    
    fprintf('First stereo..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,1)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,7)));
    [cvL2, cvR2] = tri_stereo(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 2);
    
    fprintf('Second stereo..\n');
    
    imgL = rgb2gray(uint8(EIs(:,:,:,2)));
    imgC = rgb2gray(uint8(EIs(:,:,:,4)));
    imgR = rgb2gray(uint8(EIs(:,:,:,6)));
    [cvL3, cvR3] = tri_stereo(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 3);
    
    fprintf('Third stereo..\n');
    
    [d1, dL1] = min(cvL1, [], 3);
    [d1, dR1] = min(cvR1, [], 3);
    [d1, dL2] = min(cvL2, [], 3);
    [d1, dR2] = min(cvR2, [], 3);
    [d1, dL3] = min(cvL3, [], 3);
    [d1, dR3] = min(cvR3, [], 3);
    
    i1 = dL1.*double(matte).*(dL1>1);
    i2 = dR1.*double(matte).*(dR1>1);
    i3 = dL2.*double(matte).*(dL2>1);
    i4 = dR2.*double(matte).*(dR2>1);
    i5 = dL3.*double(matte).*(dL3>1);
    i6 = dR3.*double(matte).*(dR3>1);
    
    icount = (dL1>1).*(matte)+(dL2>1).*(matte)+(dL3>1).*(matte)+(dR1>1).*(matte)+(dR2>1).*(matte)+(dR3>1).*(matte);
    
    subplot(2,3,1);
    imagesc(i1)
    subplot(2,3,2);
    imagesc(i2)
    subplot(2,3,3);
    imagesc(i3)
    subplot(2,3,4);
    imagesc(i4)
    subplot(2,3,5);
    imagesc(i5)
    subplot(2,3,6);
    imagesc(i6)
    
    figure
    i_avg = (dL1+dL2+dL3+dR1+dR2+dR3)  ./ icount.*(icount>0);
    imagesc(i_avg);
    
    
    
    cost_volume = alpha.*cost_volume_SAD+(1-alpha).*cost_volume_CENSUS;
    [d, depthAD] = min(cost_volume_SAD, [], 3);
    [d, depthCE] = min(cost_volume_CENSUS, [], 3);
    [d, depthCV] = min(cost_volume, [], 3);
    imagesc(depthAD)
    figure, imagesc(depthCE)
    figure, imagesc(depthCV)
    %{
    selpix_patch = zeros(ws*3+2, ws*3+2, size(cost_volume,3));
    cx = int32(round(selpix(1)+1));
    cy = int32(round(selpix(2)+1));
    c_img_patch = c_img(cx-hws:cx+hws, cy-hws:cy+hws);
    for l = 1:size(cost_volume,3)
        ncc_val = 0.0;
        sad_val = 0.0;
        cx = int32(round(selpix(1) - 0.866*l));
        cy = int32(round(selpix(2) - 0.5*l));
        selpix_patch(1:ws,hws+1:ws+hws,l) = big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,1);
        ncc_val = ncc_val + ncc_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,1), double(c_img_patch));
        sad_val = sad_val + sad_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,1), double(c_img_patch));
        cx = int32(round(selpix(1) - 0.866*l));
        cy = int32(round(selpix(2) + 0.5*l));
        selpix_patch(1:ws,ws+hws+2:2*ws+hws+1,l) = big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,2);
        ncc_val = ncc_val + ncc_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,2), double(c_img_patch));
        sad_val = sad_val + sad_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,2), double(c_img_patch));
        cx = int32(round(selpix(1) - 0.0));
        cy = int32(round(selpix(2) - l));
        selpix_patch(ws+2:2*ws+1,1:ws,l) = big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,3);
        ncc_val = ncc_val + ncc_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,3), double(c_img_patch));
        sad_val = sad_val + sad_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,3), double(c_img_patch));
        cx = int32(round(selpix(1)));
        cy = int32(round(selpix(2)));
        selpix_patch(ws+2:2*ws+1,ws+2:2*ws+1,l) = big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,4);
        ncc_val = ncc_val + ncc_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,4), double(c_img_patch));
        sad_val = sad_val + sad_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,4), double(c_img_patch));
        cx = int32(round(selpix(1)));
        cy = int32(round(selpix(2) + l));
        selpix_patch(ws+2:2*ws+1,2*ws+3:3*ws+2,l) = big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,5);
        ncc_val = ncc_val + ncc_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,5), double(c_img_patch));
        sad_val = sad_val + sad_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,5), double(c_img_patch));
        cx = int32(round(selpix(1) + 0.866*l));
        cy = int32(round(selpix(2) - 0.5*l));
        selpix_patch(2*ws+3:3*ws+2,hws+2:ws+hws+1,l) = big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,6);
        ncc_val = ncc_val + ncc_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,6), double(c_img_patch));
        sad_val = sad_val + sad_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,6), double(c_img_patch));
        cx = int32(round(selpix(1) + 0.866*l));
        cy = int32(round(selpix(2) + 0.5*l));
        selpix_patch(2*ws+3:3*ws+2,ws+hws+2:2*ws+hws+1,l) = big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,7);
        ncc_val = ncc_val + ncc_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,7), double(c_img_patch));
        sad_val = sad_val + sad_patch(big_matrix(cx-hws:cx+hws, cy-hws:cy+hws,7), double(c_img_patch));
        
        selpix_cost(l,1) = alpha*ncc_val + (1-alpha)*sad_val;
    end
    selpix_real_cost = cost_volume(selpix(1)+1, selpix(2)+1, :);
    lenn = 5; % size(cost_volume,3)
    %}
    %{
    figure(1)
    for i = 1:lenn
        subplot(1,lenn,i);
        imagesc(selpix_patch(:,:,i));
    end
    figure(2)
    for i = lenn+1:2*lenn
        subplot(1,lenn,i-lenn);
        imagesc(selpix_patch(:,:,i));
    end
    figure(3)
    for i = 2*lenn+1:3*lenn
        subplot(1,lenn,i-2*lenn);
        imagesc(selpix_patch(:,:,i));
    end
    figure(4)
    for i = 3*lenn+1:4*lenn
        subplot(1,lenn,i-3*lenn);
        imagesc(selpix_patch(:,:,i));
    end
    %}
    %second_row = linspace(lenn+1, lenn*2, lenn);
    %subplot(2,lenn,second_row);
    %figure()
    %plot(squeeze(selpix_real_cost));
    %hold on
    %plot(squeeze(selpix_cost));
    fprintf('Done!                       ');
    toc
end

function ncc_val = ncc_patch(patch1, patch2)

    ncc_val = correlation_coefficient(patch1, patch2);
end

function sad_val = sad_patch(patch1, patch2)

    sad_val = mean(mean(((patch1 - patch2).^2)));
end


function cc = correlation_coefficient(patch1, patch2)

    dp1 = double(patch1);
    dp2 = double(patch2);
    zmp1 = dp1 - mean(mean(dp1));
    zmp2 = dp2 - mean(mean(dp2));
    product = mean(mean(  zmp1 .* zmp2   ));
    stds = std2(dp1) * std2(dp2);
    if stds == 0
        cc = 0;
    else
        product = stds / product;
        cc =  product;
    end
end