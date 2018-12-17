function [cv, cv_sp1, cv_sp2] = depth_from_defocus_v3(fs, aif, ws, alpha, gamma_s1, gamma_s2)

    fprintf('\nComputing cost volume1 from defocus!\n')
    % compute depth from defocus by calculating difference between focal stack
    % and all-in-focus image
    tic
    num_of_imgs = size(fs,4);
    cv_sp1 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv_sp2 = zeros(size(aif,1), size(aif,2), num_of_imgs);
    cv = zeros(size(aif,1), size(aif,2), num_of_imgs);
    aif = rgb2gray(uint8(aif));
    for i=1:num_of_imgs

        fprintf('Analyzing image %d\n',i);
        colorimg = fs(:,:,:,i);
        uint8img = uint8(colorimg);
        current_fs_img = rgb2gray(uint8img);
 
        nccsadaw = compute_diff_image_aw(double(current_fs_img), double(aif), ws, double(alpha));
        aif_s1 = impyramid(aif, 'reduce');
        fs_s1 = impyramid(current_fs_img, 'reduce');
        nccsadaw_s1_tmp = compute_diff_image_aw(double(fs_s1), double(aif_s1), ws, double(alpha));
        nccsadaw_s1 = impyramid(nccsadaw_s1_tmp, 'expand');
        
        % Check size of the first scaled image (there may be a pixel
        % difference)
        slice_s1 = zeros(size(slice_s0,1), size(slice_s0,2));
        diff_s1 = [abs(size(slice_s1_tmp,1) - size(slice_s0,1)), abs(size(slice_s1_tmp,2) - size(slice_s0,2))];
        if sum(diff_s1) == 0
            % we are cool
            slice_s1 = slice_s1_tmp;
        else
            d1 = 0;
            if mod(diff_s1(1),2) ~= 0
                d1 = 1;
            end
            d2 = 0;
            if mod(diff_s1(2),2) ~= 0
                d2 = 1;
            end
            start1 = 1+floor(diff_s1(1)/2.0+d1);
            end1 = floor(diff_s1(1)/2.0);
            start2 = 1+floor(diff_s1(2)/2.0+d2);
            end2 = floor(diff_s1(2)/2.0);
            slice_s1(start1:size(slice_s1,1)-end1, start2:size(slice_s1,2)-end2) = slice_s1_tmp;
        end
        
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
        cv(:,:,i) = imguidedfilter(nccsadaw + gamma_s1.*nccsadaw_s1 + gamma_s2.*nccsadaw_s2);

    end

    fprintf('Done!                       ');
    toc
end