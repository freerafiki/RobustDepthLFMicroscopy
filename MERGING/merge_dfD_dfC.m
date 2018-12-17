function [final_depth_MRF, final_depth_FILLED] = merge_dfD_dfC(dfD_cv, dfC_cv, aif, matte, superpixels, phi, sp_size, D, WS2)

fprintf('\nMerging the two cost volumes..\n');
%% CALCULATE CONFIDENCE
fprintf('Calculating confidence..\n');
tic
% VERSION 1
[costD, depthD] = min(dfD_cv, [], 3);
[costC, depthC] = min(dfC_cv, [], 3);
[confD, confC] = IIS_computeConfidences(dfD_cv, dfC_cv, costD, costC, depthD, depthC, D, WS2);
% VERSION 2
%[dfD_conf, dfD_depth] = est_conf(dfD_cv); 
%[dfC_conf, dfC_depth] = est_conf(dfC_cv);

% SHOW CONFIDENCES
%subplot(1,2,1);
%imagesc(dfD_conf);
%subplot(1,2,2)
%imagesc(dfC_conf)

%% CONFIDENCE ASSESSMENT
sum_confD = sum(sum(confD.*double(matte))); % sum of values wihtout masked part
num_of_valsD = sum(sum(matte));
mean_confD = sum_confD / num_of_valsD;
sum_confC = sum(sum(confC.*double(matte))); % sum of values wihtout masked part
num_of_valsC = sum(sum(matte));
mean_confC = sum_confC / num_of_valsC;
%std_confD = std(std(confD));
%std_confC = std(std(confC));


%% PICK PHI
if phi == -1
    % choose phi
    phi = mean_confD / (mean_confD + mean_confC);
end
% mean confidence
% = mean_confD^2 + mean_confC^2 / (mean_confD + mean_confC)
mean_conf = phi*mean_confD + (1-phi)*mean_confC; 
confMIX = (phi*confD + (1-phi)*confC).*double(matte);


%% USING SUPERPIXELS/GRID AND A OCCWEIGHT SIMILAR IDEA
% sum the two contributions and create an "occlusion" or "outliers" map
% that will be filled using the reliable values
fprintf('First method..\n');
MIN_RATIO_CONF = 0.9 / mean_conf; % so that first round is checking conf > 0.9!
STEP_REDUCTION = 0.1;
MIN_RATIO_FILLED = 0.45;
AREA_NO_MATTED = sum(sum(matte>0));
[~, final_depth_TOFILL] = min(phi.*dfD_cv + (1-phi).*dfC_cv, [], 3);
final_depth_TOFILL = final_depth_TOFILL .*double(matte);
final_depth_FILLED = final_depth_TOFILL; %zeros(size(final_depth_TOFILL,1), size(final_depth_TOFILL,2));

%{
%% SUPERPIXELS VERSION 1 

fprintf('Finding outliers..\n');
THRESH_RATIO_CONF = 1.0;
THRESH_RATIO_FILLED = 0.2;
CONF_RATIO_STEP = 0.1;
outliers_finished = false;
while (outliers_finished == false)
    [row, col, v] = find(confMIX < THRESH_RATIO_CONF*mean_conf & confMIX > 0);
    ratio = sum(v) / sum(sum(matte));
    if ratio < THRESH_RATIO_FILLED
        outliers_finished = true;
    else
        THRESH_RATIO_CONF = THRESH_RATIO_CONF - CONF_RATIO_STEP;
    end
end

fprintf('Filling the values..\n');
% half window size (hws)
hws = floor(WS2/2.0);
g1 = WS2 / 2;
g2 = mean(max(max(aif))) / 10;
g3 = 10;
for j = 1: size(row,1)
    % row and col are the indices of the outliers - pixels to be replaced
    x = row(j);
    y = col(j);
    if v > 0
        % calculate likelihood of the neighbours
        max_likelihood = 0;
        depth_at_max_likelihood_point = 0;
        for h=-hws:+hws
           for w=-hws:+hws
               spatial_dist_factor = sqrt(h^2+w^2);
               if x+h<size(aif,1) | y+w<size(aif,2)
                   color_dist_factor = g2;
                   conf_factor = g3;
               else
                   color_dist_factor = mean(abs(aif(x+h, y+w)-aif(x,y)));
                   conf_factor = confMIX(x+h,y+w);
               end
               
               likelihood = exp((-(spatial_dist_factor/g1)-(color_dist_factor/g2)-(conf_factor/g3)));
               if likelihood > max_likelihood
                   max_likelihood = likelihood;
                   depth_at_max_likelihood_point = final_depth_TOFILL(x,y);
               end
           end
        end
        final_depth_FILLED(x,y) = depth_at_max_likelihood_point;
        subplot(121)
        imagesc(final_depth_TOFILL(x-5:x+5, y-5:y+5))
        subplot(122)
        imagesc(final_depth_FILLED(x-5:x+5, y-5:y+5))
    end
end

imagesc(final_depth_TOFILL), figure, imagesc(final_depth_FILLED)

imagesc((confMIX < mean_conf) .* (confMIX > 0) .* matte);
sum(v)
sum(sum(matte))
sp_max_val = max(max(superpixels));
for i = 1: sp_max_val
    % get pixels
    indices = superpixels == i;
    if sum(sum(indices.*matte)) > MIN_NUM_PIXELS
        curRegionDFD = depthD.*indices;
        curRegionDFC = depthC.*indices;
        curRegionConfDFD = confD.*indices;
        curRegionConfDFC = confC.*indices;
        confValuesDFD = curRegionConfDFC(superpixels==i);
        confValuesDFC = curRegionConfDFC(superpixels==i);
        validConfValuesDFD = confValuesDFD(confValuesDFD>0.00);
        validConfValuesDFC = confValuesDFC(confValuesDFC>0.00);
        muDFD = mean((validConfValuesDFD));
        muDFC = mean((validConfValuesDFC));
        % pick the bad pixels
        curRegionThresh = phi*muDFD + (1-phi)*muDFC;
        curRegionConfMIX = phi.*curRegionConfDFD+(1-phi).*curRegionConfDFC;
        outliers = curRegionConfMIX < curRegionThresh;
    end
end

%% SUPERPIXELS LATER VERSION INCLUDING NEIGHBOURS - NOT READY YET

sp_max_val = max(max(superpixels));
% loop for conversion
first_step = false;
good_sp = zeros(size(superpixels,1), size(superpixels,2));
while (first_step == false);
    % mark good superpixels;
    for i = 1: sp_max_val
        % get pixels
        indices = superpixels == i;
        if sum(sum(indices.*matte)) > MIN_NUM_PIXELS
            curRegionConfDFD = confD.*indices;
            curRegionConfDFC = confC.*indices;
            confValuesDFD = curRegionConfDFD(superpixels==i);
            confValuesDFC = curRegionConfDFC(superpixels==i);
            validConfValuesDFD = confValuesDFD(confValuesDFD>0.00);
            validConfValuesDFC = confValuesDFC(confValuesDFC>0.00);
            muDFD = mean((validConfValuesDFD));
            muDFC = mean((validConfValuesDFC));
            % select good pixels
            if (phi)*muDFD + (1-phi)*muDFC > MIN_RATIO_CONF*mean_conf
                good_sp = good_sp + indices.*matte;
            end
        end
    end
    if sum(sum(good_sp)) > AREA_NO_MATTED * MIN_RATIO_FILLED;
        first_step = true;
    else
        MIN_RATIO_CONF = MIN_RATIO_CONF - STEP_REDUCTION;
    end
end              
%}

%% USING MRF from TAO
% Markov random field approach to minimize energy that has a term for each
% cost volume (correspondences and defocus) - parameters changed to obtain
% less smoothed results
addpath(genpath('TAO_COMBINING2013'));
fprintf('Second Method..\n');
%%% Regularize                          --------------
WS_PENALTY_W1           = 0.3                                             ;
WS_PENALTY_W2           = 0.1                                             ;
lambda_flat             = 2                                               ;
lambda_smooth           = 2                                               ;
ROBUSTIFY_SMOOTHNESS    = 1                                               ;
gradient_thres          = 1.0                                             ;
SOFTEN_EPSILON          = 1.0                                             ;
CONVERGE_FRACTION       = 1;  
% GATHER PARAMETERS
parameters       = struct('WS_PENALTY_W1',WS_PENALTY_W1,...
                             'WS_PENALTY_W2',WS_PENALTY_W2,...
                             'lambda_flat',lambda_flat,...
                             'lambda_smooth',lambda_smooth,...
                             'ROBUSTIFY_SMOOTHNESS',ROBUSTIFY_SMOOTHNESS,...
                             'gradient_thres',gradient_thres,...
                             'SOFTEN_EPSILON',SOFTEN_EPSILON,...
                             'CONVERGE_FRACTION',CONVERGE_FRACTION...
                             ) ;
final_depth_MRF        = DEPTH_MRF(depthD, depthC,...
                                confD, confC,...
                                aif,parameters);
                            
final_depth_MRF = final_depth_MRF .* double(matte);
fprintf('Done!                       ');
toc