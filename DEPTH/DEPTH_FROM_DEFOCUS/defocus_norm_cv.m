function [cv, cv_sp1, cv_sp2] = defocus_norm_cv(fs, aif, ws, alpha, gamma1, gamma2, gamma3, priors)

    fprintf('\nComputing cost volume1 from defocus!\n')
    % compute depth from defocus by calculating difference between focal stack
    % and all-in-focus image
    tic
    num_of_imgs = size(fs,4);
    cv_sp1 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv_sp2 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv = zeros(size(aif,1), size(aif,2), num_of_imgs);
    
    % THIS IS PICTURE DEPENDENT
    % 1 = FiMic
    % 2 = Lytro
    PICTURE_TYPE = 1;
    if PICTURE_TYPE == 1
        % FOR FiMic needs uint8 conversion
        gaif = rgb2gray(uint8(aif)); % FIMIC VERSION
    elseif PICTURE_TYPE == 2
        % FOR LYTRO NOT NEEDED
        gaif = rgb2gray(aif);
    end
    %cd DEPTH/CPP/
    %mex calc_ncc.cpp;
    %mex calc_sad.cpp;
    %mex calc_diff.cpp;
    %mex compute_diff_image.cpp
    %mex compute_diff_image_aw.cpp
    %cd ../..
    %gamma_s1 = 0.6;
    %gamma_s2 = 0.3;
    %hws = floor(ws/2);
    %selpix = [150, 364];
    %selpix_patch = zeros(ws, ws, num_of_imgs);
    %selpix_cost = zeros(num_of_imgs,1);
    %selpix_real_cost = zeros(num_of_imgs,1);
    
    %% PREPARE PRIORS MAPS
    thresh = multithresh(priors,2);
    maps = imquantize(priors,thresh);
    high_freqs = maps == 3;
    med_freqs = maps == 2;
    low_freqs = maps == 1;
    
    %% PRECOMPUTE GUIDED FILTER
    r = 5;
    eps = 0.001;
    guidedfilter_color_precompute(aif, r, eps);
    
    for i=1:num_of_imgs

        fprintf('Analyzing image %d\n',i);
        colorimg = fs(:,:,:,i);
        if PICTURE_TYPE == 1
            % FOR FiMic needs uint8 conversion
            % FIMIC VERSION
            uint8img = uint8(colorimg);
            current_fs_img = rgb2gray(uint8img);
        elseif PICTURE_TYPE == 2
            % FOR LYTRO NOT NEEDED
            current_fs_img = rgb2gray(colorimg);
        end
        
        %% VERSION 2 - calculate difference normally and then choose maximum using superpixels
        %tic
        %ncc = calc_ncc(double(current_fs_img), double(aif), ws);
        %sad = calc_sad(double(current_fs_img), double(aif), ws);
        %nccsad = compute_diff_image(double(current_fs_img), double(aif), ws, double(alpha));
        nccsadaw = compute_diff_image_aw(double(current_fs_img), double(gaif), ws, double(alpha));
        aif_s1 = impyramid(gaif, 'reduce');
        fs_s1 = impyramid(current_fs_img, 'reduce');
        nccsadaw_s1_tmp = compute_diff_image_aw(double(fs_s1), double(aif_s1), ws, double(alpha));
        nccsadaw_s1 = impyramid(nccsadaw_s1_tmp, 'expand');
        
        fs_s2 = impyramid(fs_s1, 'reduce');
        aif_s2 = impyramid(aif_s1, 'reduce');
        nccsadaw_s2_diff = compute_diff_image_aw(double(fs_s2), double(aif_s2), ws, double(alpha));
        nccsadaw_s2_tmp = impyramid(impyramid(nccsadaw_s2_diff, 'expand'), 'expand');
        if size(nccsadaw_s2_tmp,1) < size(nccsadaw_s1,1)
            nccsadaw_s2 = zeros(size(nccsadaw_s1,2), size(nccsadaw_s1,2));
            nccsadaw_s2(2:size(nccsadaw_s2,1)-1, 2:size(nccsadaw_s2,2)-1) = nccsadaw_s2_tmp;
        else
            nccsadaw_s2 = nccsadaw_s2_tmp;
        end
        
        cur_image = high_freqs .* (nccsadaw .* gamma1 + nccsadaw_s1 .* gamma2 + nccsadaw_s2 .* gamma3) + ...
            med_freqs .* (nccsadaw .* gamma3 + nccsadaw_s1 .* gamma1 + nccsadaw_s2 .* gamma2) + ...
            low_freqs .* (nccsadaw .* gamma3 + nccsadaw_s1 .* gamma2 + nccsadaw_s2 .* gamma1);
      	
        %cur_image2 = (nccsadaw + gamma_s1.*nccsadaw_s1 + gamma_s2.*nccsadaw_s2) / (1 + gamma_s1 + gamma_s2);
        cv(:,:,i) = min(1, max(0, guidedfilter_color_runfilter(cur_image)));
        %toc

    end
    
    
    fprintf('Done!                       ');
    toc
end
