function cv_proc = invert_volume(cv_ini)

%% INVERT THE VOLUME
% the disparity from correspondences and defocus have opposite order
% so that a minimum of x in defocus in N-x in correpsondences
% where N is the max disparity
cv_proc = zeros(size(cv_ini,1), size(cv_ini,2), size(cv_ini,3));
for j = 1:size(cv_ini,3)
    cv_proc(:,:,j) = cv_ini(:,:,size(cv_ini,3)-j+1);
end
    
