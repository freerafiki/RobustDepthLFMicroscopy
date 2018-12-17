function [cv, cv_sp1, cv_sp2, cv_sp3] = depth_from_defocus(fs, aif, ws, alpha)

    fprintf('\nComputing cost volume1 from defocus!\n')
    % compute depth from defocus by calculating difference between focal stack
    % and all-in-focus image
    tic
    num_of_imgs = size(fs,4);
    cv_sp1 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv_sp2 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv_sp3 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv2 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    aif = rgb2gray(uint8(aif));
    cd DEPTH/CPP/
    mex calc_ncc.cpp;
    mex calc_sad.cpp;
    %mex calc_diff.cpp;
    mex compute_diff_image.cpp
    cd ../..
    hws = floor(ws/2);
    selpix = [150, 364];
    selpix_patch = zeros(ws, ws, num_of_imgs);
    selpix_cost = zeros(num_of_imgs,1);
    selpix_real_cost = zeros(num_of_imgs,1);
    for i=1:num_of_imgs

        fprintf('Analyzing image %d\n',i);
        colorimg = fs(:,:,:,i);
        uint8img = uint8(colorimg);
        current_fs_img = rgb2gray(uint8img);
        
        %% VERSION 1 - calculate difference using superpixels
        
        %tic
        %[sp_big, sp_small, diff] = calc_diff(double(current_fs_img), double(aif), double(superpixels1), double(superpixels2), double(matte), ws, alpha);
        %cv_sp1(:,:,i) = imguidedfilter(sp_big);
        %cv_sp2(:,:,i) = imguidedfilter(sp_small);
        %cv(:,:,i) = imguidedfilter(diff,'NeighborhoodSize',[7 7]);
        %toc
        
        
        %% VERSION 2 - calculate difference normally and then choose maximum using superpixels
        %tic
        ncc = calc_ncc(double(current_fs_img), double(aif), ws);
        sad = calc_sad(double(current_fs_img), double(aif), ws);
        nccsad = compute_diff_image(double(current_fs_img), double(aif), ws, double(alpha));
        cur_diff_image = ncc*alpha + (1-alpha)*sad;
        cv(:,:,i) = imguidedfilter(nccsad); %,'NeighborhoodSize',[7 7]);
        cv_sp1(:,:,i) = imguidedfilter(cur_diff_image);
        cv_sp2(:,:,i) = imguidedfilter(ncc);
        cv_sp3(:,:,i) = imguidedfilter(sad);
        %toc
        
        %% DEBUG VISUALIZATION
        %cx = selpix(1);
        %cy = selpix(2);
        %selpix_patch(:,:,i) = current_fs_img(cx-hws:cx+hws, cy-hws:cy+hws);
        %aif_patch = aif(cx-hws:cx+hws, cy-hws:cy+hws);
        %ncc_val = ncc_patch(selpix_patch(:,:,i), double(aif_patch));
        %sad_val = sad_patch(selpix_patch(:,:,i), double(aif_patch));
        %selpix_cost(i,1) = alpha*ncc_val + (1-alpha)*sad_val;
        %selpix_real_cost(i,1) = cur_diff_image(cx, cy);
        
        %cv(:,:,i) = cur_diff_image;
        %{
        min_sp = min(min(superpixels));
        max_sp = max(max(superpixels));

        for cc = min_sp:max_sp
            fprintf('cc %d\n', cc);
            sp_img = (superpixels==cc) .* logical(matte);
            if sum(sum(sp_img)) > 0
                avg_val = sum(sum(cv(:,:,i).*sp_img))/sum(sum(sp_img));
                cv_sp1 = cv_sp1(:,:,i) + avg_val.*sp_img;
            end
        end
        %}
        %ncc2 = calc_ncc(double(current_fs_img), double(aif));
        %ncc = ncc / max(max(ncc));
        %ncc = ncc/2;
        %subplot(2,2,1)
        %imshow(current_fs_img)
        %subplot(2,2,2)
        %imshow(aif)
        %subplot(2,2,3)
        %imagesc(ncc.*double(matte), [0,30])
        %subplot(2,2,4)
        %imagesc(sad.*double(matte), [0,1])
        
        
        %sad_diff(:,:,i) = sad;
    end
    
    
    
    
    
    %% DEBUG PLOT
    %{
    for i = 1:num_of_imgs
        subplot(3,num_of_imgs,i);
        imagesc(selpix_patch(:,:,i));
        subplot(3,num_of_imgs,num_of_imgs+i);
        imagesc(aif_patch);
    end
    third_row = linspace(num_of_imgs*2+1, num_of_imgs*3, num_of_imgs);
    subplot(3,num_of_imgs,third_row)
    plot(selpix_cost);
    figure(); % hold on
    plot(selpix_real_cost);
    %}
    
    %% VERSION 1 - pick minimum
    %[values, disp] = min(cv_sp1,[],3);
    %[values2, disp2] = min(cv_sp2,[],3);
    %[values3, disp3] = min(cv,[],3);
    
    %% VERSION 2 - Pick minimum using superpixels
    %cd DEPTH/CPP/
    %mex extract_argmin.cpp
    %cd ../..
    %[depth_sp_big, depth_sp_small, depth] = extract_argmin(cv2, double(superpixels1), double(superpixels2), double(matte));
    
    %map to real depth value
    %{
    fl = zeros(size(fs,1));
    for i=1:size(fs,1)
        fl(i) = fs{i,2};
    end
    mex map_fl_values.cpp
    real_depth = map_fl_values(depth, fl);
    
    %depth = depth_vals(disp);
    %depth2 = depth_vals(disp2);
    subplot(1,3,1)
    imagesc(depth, [0 20])
    subplot(1,3,2)
    imagesc(depth_sp_small, [0 20])
    subplot(1,3,3)
    imagesc(depth_sp_big, [0 20])
    hold on
    colormap(jet)
    %}
    
    
    fprintf('Done!                       ');
    toc
end

function ncc_val = ncc_patch(patch1, patch2)

    ncc_val = correlation_coefficient(patch1, patch2);
end

function sad_val = sad_patch(patch1, patch2)

    sad_val = mean(mean(((patch1 - patch2).^2)));
end

function ncc = calculate_ncc(img, ref)

    ncc = zeros(size(img,1), size(img,2));
    d = 3;
    for k=d+1:(size(img,1)-d) 
        for l=d+1:(size(img,2)-d)
            if sqrt((k-size(img,1)/2)^2 + (l-size(img,2)/2)^2) < 400
                %fprintf('Calculating pixel (%d, %d)\n',k,l);
                ncc(k, l) = correlation_coefficient(ref(k-d:k+d, l-d:l+d), img(k-d:k+d, l-d:l+d));
            end
        end
    end
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
        product = product/ stds;
        cc =  product;
    end
end