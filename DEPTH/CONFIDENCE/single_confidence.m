function conf = single_confidence(cv, img, depth, which_conf, matte)

%% COMPUTE CONFIDENCE
% check paper "Hu et al. - A Quantitative Evaluation of Confidence Measures
% for Stereo Vision, PAMI 2012" for more info on measures
% we use only measures that look at the cost curve (group 1 in the paper) 
% to maintain coherence between defocus and correspondence cost cubes 
% (stereo would have the left and right costs, but defocus wouldn't)


conf = ones(size(img,1), size(img,2));

if which_conf == 0
    %% MSM
    conf = -cv(:,:,depth);
elseif which_conf == 1
    %% PKR  --> modified because our depth is computed with MGM and not WTA
    % idea here is ratio between chosen value and second smallest minimum
    % could be optimized using matrices and not "for" loops
    % conf_D(:,:) = min(1, unique(findPeaks(cv_D(:,:,:)))(2) / cv_D(:,:,depth_D(:,:)) );
    % ???
    for i = 1:size(img,1)
        for j = 1:size(img,2)
            if matte(i,j) == 1
                peaks = findPeaks(cv(i,j));
                sortedPeaks = unique(peaks, 'first');
                conf(i,j) = min(1, sortedPeaks(2) / cv(i,j,depth(i,j)) );
            else
                conf(i,j) = 0;
            end
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
            if matte(i,j) == 1
                values = unique(cv(i,j,:));
                if numel(values) == 1
                    conf(i,j) = 1;
                else
                    conf(i,j) = min(10, max(1, values(2) / cv(i,j,depth(i,j)) ) );
                end
            else
                conf(i,j) = 0;
            end
        end
    end
    % NORMALIZE
    conf = conf / max(max(conf));
elseif which_conf == 3
    %% AML - Attainable Maximum Likelihood
    % conf_D(:,:) = 1 / pow(sum((cv_D(:,:,:) - depth_D,3),2) ..
    sigma_AML = 1;
    for i = 1:size(img,1)
        for j = 1:size(img,2)
            if matte(i,j) == 1
                denomD = sum( exp(-pow(cv(i,j,:)-depth(i,j),2)/2 * pow(sigma_AML,2)) );
                conf(i,j) = 1 / denomD;
            else
                conf(i,j) = 0;
            end
        end
    end
end

%% accounts for priors
use_priors = false;
if use_priors
    factor = 0.8;
    conf = (factor).*conf + (1-factor).*priorsmap;
end

%% apply matte
conf = conf .* double(matte);