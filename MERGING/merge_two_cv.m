function  final_cost_volume = merge_two_cv(cv1, cv2, matte, D, WS2)

% VERSION 1
[costD, depthD] = min(cv1, [], 3);
[costC, depthC] = min(cv2, [], 3);
[confD, confC] = IIS_computeConfidences(cv1, cv2, costD, costC, depthD, depthC, D, WS2);
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
phi = mean_confD / (mean_confD + mean_confC);
% mean confidence
% = mean_confD^2 + mean_confC^2 / (mean_confD + mean_confC)
mean_conf = phi*mean_confD + (1-phi)*mean_confC; 
confMIX = (phi*confD + (1-phi)*confC).*double(matte);
final_cost_volume = phi.*cv1 + (1-phi).*cv2;
