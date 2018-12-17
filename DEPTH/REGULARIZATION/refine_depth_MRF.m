function final_depth = refine_depth_MRF(depth_D, depth_C, conf_D, conf_C, img, matte)

%% USING MRF from TAO
% Markov random field approach to minimize energy that has a term for each
% cost volume (correspondences and defocus) - parameters changed to obtain
% less smoothed results
tic

%addpath(genpath('TAO_COMBINING2013'));
fprintf('Regularization..\n');
%%% Regularize                          --------------
WS_PENALTY_W1           = 10                                             ;
WS_PENALTY_W2           = 20                                             ;
lambda_flat             = 2                                               ;
lambda_smooth           = 1                                               ;
ROBUSTIFY_SMOOTHNESS    = 1                                               ;
gradient_thres          = 1.0                                             ;
SOFTEN_EPSILON          = 1.0                                             ;
CONVERGE_FRACTION       = 0.5;  
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
final_depth_MRF        = DEPTH_MRF(depth_D, depth_C,...
                                conf_D, conf_C,...
                                img,parameters);
                            
final_depth = final_depth_MRF .* double(matte);
fprintf('Done!                       ');
toc