function [cvMS, cvMSGF, cvMSBF] = depth_from_correspondences_v3(EIs, c_img, ws, alpha, d_min, d_max, matte, gamma_s1, gamma_s2, index)

    fprintf('\nComputing depth from correspondences!\n')
    % computes disparity from correspondences using multi view approach
    tic

    %% CALCULATE STEREO USING THE THREE EPIPOLAR LINES
    
    %% FIRST EPIPOLAR LINE
    fprintf('First Epipolar Line..\n');
    
    % choose case for image
    
    % If it's the central use the tri_stereo (use three images)
    
    if index == 4
        
        imgL = rgb2gray(uint8(EIs(:,:,:,3)));
        imgC = rgb2gray(uint8(EIs(:,:,:,4)));
        imgR = rgb2gray(uint8(EIs(:,:,:,5)));
        [cvL1, cvR1] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 1);
        imgL_s1 = impyramid(imgL, 'reduce');
        imgL_s2 = impyramid(imgL_s1, 'reduce');
        imgC_s1 = impyramid(imgC, 'reduce');
        imgC_s2 = impyramid(imgC_s1, 'reduce');
        imgR_s1 = impyramid(imgR, 'reduce');
        imgR_s2 = impyramid(imgR_s1, 'reduce');
        [cvL1_s1, cvR1_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 1);
        [cvL1_s2, cvR1_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 1);
    
    else
        % for all other images at the border pick only two 
        % this is only first epipolar line
        % imgL is always the reference image
        % direction controls the dmin and dmax (use positive or negative
        % disparities)
        if index == 1
        
            imgL = rgb2gray(uint8(EIs(:,:,:,1)));
            imgR = rgb2gray(uint8(EIs(:,:,:,2)));
            direction = 1;
        
        elseif index == 2
            
            imgL = rgb2gray(uint8(EIs(:,:,:,2)));
            imgR = rgb2gray(uint8(EIs(:,:,:,1)));
            direction = -1;

        elseif index == 3
                      
            imgL = rgb2gray(uint8(EIs(:,:,:,3)));
            imgR = rgb2gray(uint8(EIs(:,:,:,4)));
            direction = 1;

        elseif index == 5
            
            imgL = rgb2gray(uint8(EIs(:,:,:,5)));
            imgR = rgb2gray(uint8(EIs(:,:,:,4)));
            direction = -1;
            
        elseif index == 6
            
            imgL = rgb2gray(uint8(EIs(:,:,:,6)));
            imgR = rgb2gray(uint8(EIs(:,:,:,7)));
            direction = 1;
        elseif index == 7
            
            imgL = rgb2gray(uint8(EIs(:,:,:,7)));
            imgR = rgb2gray(uint8(EIs(:,:,:,6)));
            direction = -1;
            
        end   
        
        cvL1 = bin_stereo_v3(double(imgL), double(imgR), (d_min), (d_max), (ws), (alpha), 1, direction);
        cvR1 = cvL1;
        imgL_s1 = impyramid(imgL, 'reduce');
        imgL_s2 = impyramid(imgL_s1, 'reduce');
        imgR_s1 = impyramid(imgR, 'reduce');
        imgR_s2 = impyramid(imgR_s1, 'reduce');
        cvL1_s1 = bin_stereo_v3(double(imgL_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 1, direction);
        cvL1_s2 = bin_stereo_v3(double(imgL_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 1, direction);
        cvR1_s1 = cvL1_s1;
        cvR1_s2 = cvL1_s2;
        
    end
    
    %% SECOND EPIPOLAR LINE
    fprintf('Second Epipolar Line..\n');
    
    if index == 4
        
        imgL = rgb2gray(uint8(EIs(:,:,:,1)));
        imgC = rgb2gray(uint8(EIs(:,:,:,4)));
        imgR = rgb2gray(uint8(EIs(:,:,:,7)));
        [cvL2, cvR2] = tri_stereo(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 2);
        imgL_s1 = impyramid(imgL, 'reduce');
        imgL_s2 = impyramid(imgL_s1, 'reduce');
        imgC_s1 = impyramid(imgC, 'reduce');
        imgC_s2 = impyramid(imgC_s1, 'reduce');
        imgR_s1 = impyramid(imgR, 'reduce');
        imgR_s2 = impyramid(imgR_s1, 'reduce');
        [cvL2_s1, cvR2_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 2);
        [cvL2_s2, cvR2_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 2);

    else
        % for all other images at the border pick only two 
        % this is only first epipolar line
        % imgL is always the reference image
        % direction controls the dmin and dmax (use positive or negative
        % disparities)
        if index == 1
        
            imgL = rgb2gray(uint8(EIs(:,:,:,1)));
            imgR = rgb2gray(uint8(EIs(:,:,:,4)));
            direction = 1;
        
        elseif index == 2
            
            imgL = rgb2gray(uint8(EIs(:,:,:,2)));
            imgR = rgb2gray(uint8(EIs(:,:,:,5)));
            direction = 1;

        elseif index == 3
                      
            imgL = rgb2gray(uint8(EIs(:,:,:,3)));
            imgR = rgb2gray(uint8(EIs(:,:,:,6)));
            direction = 1;

        elseif index == 5
            
            imgL = rgb2gray(uint8(EIs(:,:,:,5)));
            imgR = rgb2gray(uint8(EIs(:,:,:,2)));
            direction = -1;
            
        elseif index == 6
            
            imgL = rgb2gray(uint8(EIs(:,:,:,6)));
            imgR = rgb2gray(uint8(EIs(:,:,:,3)));
            direction = -1;
            
        elseif index == 7
            
            imgL = rgb2gray(uint8(EIs(:,:,:,7)));
            imgR = rgb2gray(uint8(EIs(:,:,:,4)));
            direction = -1;
            
        end 
        
        cvL2 = bin_stereo_v3(double(imgL), double(imgR), (d_min), (d_max), (ws), (alpha), 1, direction);
        cvR2 = cvL2;
        imgL_s1 = impyramid(imgL, 'reduce');
        imgL_s2 = impyramid(imgL_s1, 'reduce');
        imgR_s1 = impyramid(imgR, 'reduce');
        imgR_s2 = impyramid(imgR_s1, 'reduce');
        cvL2_s1 = bin_stereo_v3(double(imgL_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 2, direction);
        cvL2_s2 = bin_stereo_v3(double(imgL_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 2, direction);
        cvR2_s1 = cvL2_s1;
        cvR2_s2 = cvL2_s2;
        
    end
    
    %% THIRD EPIPOLAR LINE
    fprintf('Third Epipolar Line..\n');
    
    if index == 4
        
        imgL = rgb2gray(uint8(EIs(:,:,:,2)));
        imgC = rgb2gray(uint8(EIs(:,:,:,4)));
        imgR = rgb2gray(uint8(EIs(:,:,:,6)));
        [cvL3, cvR3] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), (d_min), (d_max), (ws), (alpha), 3);
        imgL_s1 = impyramid(imgL, 'reduce');
        imgL_s2 = impyramid(imgL_s1, 'reduce');
        imgC_s1 = impyramid(imgC, 'reduce');
        imgC_s2 = impyramid(imgC_s1, 'reduce');
        imgR_s1 = impyramid(imgR, 'reduce');
        imgR_s2 = impyramid(imgR_s1, 'reduce');
        [cvL3_s1, cvR3_s1] = tri_stereo_v2(double(imgL_s1), double(imgC_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 3);
        [cvL3_s2, cvR3_s2] = tri_stereo_v2(double(imgL_s2), double(imgC_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 3);
    
    else
        % for all other images at the border pick only two 
        % this is only first epipolar line
        % imgL is always the reference image
        % direction controls the dmin and dmax (use positive or negative
        % disparities)
        if index == 1
        
            imgL = rgb2gray(uint8(EIs(:,:,:,1)));
            imgR = rgb2gray(uint8(EIs(:,:,:,3)));
            direction = -1;
        
        elseif index == 2
            
            imgL = rgb2gray(uint8(EIs(:,:,:,2)));
            imgR = rgb2gray(uint8(EIs(:,:,:,4)));
            direction = -1;

        elseif index == 3
                      
            imgL = rgb2gray(uint8(EIs(:,:,:,3)));
            imgR = rgb2gray(uint8(EIs(:,:,:,1)));
            direction = 1;

        elseif index == 5
            
            imgL = rgb2gray(uint8(EIs(:,:,:,5)));
            imgR = rgb2gray(uint8(EIs(:,:,:,7)));
            direction = -1;
            
        elseif index == 6
            
            imgL = rgb2gray(uint8(EIs(:,:,:,6)));
            imgR = rgb2gray(uint8(EIs(:,:,:,4)));
            direction = 1;
            
        elseif index == 7
            
            imgL = rgb2gray(uint8(EIs(:,:,:,7)));
            imgR = rgb2gray(uint8(EIs(:,:,:,5)));
            direction = +1;
            
        end 
        
        cvL3 = bin_stereo_v3(double(imgL), double(imgR), (d_min), (d_max), (ws), (alpha), 1, direction);
        cvR3 = cvL3;
        imgL_s1 = impyramid(imgL, 'reduce');
        imgL_s2 = impyramid(imgL_s1, 'reduce');
        imgR_s1 = impyramid(imgR, 'reduce');
        imgR_s2 = impyramid(imgR_s1, 'reduce');
        cvL3_s1 = bin_stereo_v3(double(imgL_s1), double(imgR_s1), (d_min/2), (round(d_max/2)), (round(ceil(ws/2))), (alpha), 2, direction);
        cvL3_s2 = bin_stereo_v3(double(imgL_s2), double(imgR_s2), (d_min/4), (round(d_max/4)), (round(ceil(ws/4))), (alpha), 2, direction);
        cvR3_s1 = cvL3_s1;
        cvR3_s2 = cvL3_s2;
        
    end
    
    %% MULTI - SCALE COST AGGREGATION
    fprintf('Summing up contributions using multi-scale..\n');

    cvMS = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
    cvMSGF = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
    cvMSBF = zeros(size(cvL1, 1), size(cvL1, 2), size(cvL1, 3));
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
        cvMSBF(:,:,m) = jointBF(cvMS(:,:,m), EIs(:,:,:,index), 5);
    end
    
    
    fprintf('Done!                       ');
    toc
end
