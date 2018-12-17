function [conf_D, conf_C] = compute_confidence(depth_D, cv_D, depth_C, cv_C, img, matte, disp_range, window_size, priorsmap, use_priors, which_conf)


%% COMPUTE CONFIDENCE
% check paper "Hu et al. - A Quantitative Evaluation of Confidence Measures
% for Stereo Vision, PAMI 2012" for more info on measures
% we use only measures that look at the cost curve (group 1 in the paper) 
% to maintain coherence between defocus and correspondence cost cubes 
% (stereo would have the left and right costs, but defocus wouldn't)

if min(min(depth_D)) == 0
    depth_D = depth_D +1;
end
if min(min(depth_C)) == 0
    depth_C = depth_C +1;
end

conf_D = ones(size(img,1), size(img,2));
conf_C = ones(size(img,1), size(img,2));
if which_conf == 0
    %% MSM
    conf_D = -cv_D(:,:,depth_D);
    conf_C = -cv_C(:,:,depth_C);
elseif which_conf == 1
    %% PKR  --> modified because our depth is computed with MGM and not WTA
    % idea here is ratio between chosen value and second smallest minimum
    % could be optimized using matrices and not "for" loops
    % conf_D(:,:) = min(1, unique(findPeaks(cv_D(:,:,:)))(2) / cv_D(:,:,depth_D(:,:)) );
    % ???
    for i = 1:size(img,1)
        for j = 1:size(img,2)
            peaksD = findPeaks(cv_D(i,j));
            sortedPeaksD = unique(peaksD, 'first');
            conf_D(i,j) = min(1, sortedPeaksD(2) / cv_D(i,j,depth_D(i,j)) );
            peaksC = findPeaks(cv_D(i,j));
            sortedPeaksC = unique(peaksC, 'first');
            conf_C(i,j) = min(1, sortedPeaksC(2) / cv_C(i,j,depth_C(i,j)) );
        end
    end
elseif which_conf == 2
    %% PKRN  --> modified because our depth is computed with MGM and not WTA
    % idea here is ratio between chosen value and second smallest value
    % (may be not a minimum)
    % could be optimized using matrices and not "for" loops
    % conf_D(:,:) = min(1, unique(cv_D(:,:,:))(2) / cv_D(:,:,depth_D(:,:)) );
    % ???
    for i = 1:size(img,1)
        for j = 1:size(img,2)
            valuesD = unique(cv_D(i,j,:));
            if numel(valuesD) == 1
                conf_D(i,j) = 1;
            else
                conf_D(i,j) = min(10, max(1, valuesD(2) / cv_D(i,j,depth_D(i,j)) ) );
            end
            valuesC = unique(cv_C(i,j,:));
            if numel(valuesC) == 1
                conf_C(i,j) = 1;
            else
                conf_C(i,j) = max(1, valuesC(2) / cv_C(i,j,depth_C(i,j)) );
            end
        end
    end
    % NORMALIZE
    conf_D = conf_D / max(max(conf_D));
    conf_C = conf_C / max(max(conf_C));
elseif which_conf == 3
    %% AML - Attainable Maximum Likelihood
    % conf_D(:,:) = 1 / pow(sum((cv_D(:,:,:) - depth_D,3),2) ..
    for i = 1:size(img,1)
        for j = 1:size(img,2)
            denomD = sum( exp(-pow(cv_D(i,j,:)-depth_D(i,j),2)/2 * pow(sigma_AML,2)) );
            conf_D(i,j) = 1 / denomD;
            denomC = sum( exp(-pow(cv_C(i,j,:)-depth_C(i,j),2)/2 * pow(sigma_AML,2)) );
            conf_C(i,j) = 1 / denomC;
        end
    end
end

%% accounts for priors
if use_priors
    factor = 0.8;
    conf_D = (factor).*conf_D + (1-factor).*priorsmap;
    conf_C = (factor).*conf_C + (1-factor).*priorsmap;
end

%% apply matte
conf_D = conf_D .* double(matte);
conf_C = conf_C .* double(matte);
