% TEST STEREO

%addpath(genpath('DEPTH'));
%addpath(genpath('DENOISING'));
imgL = imread('/Users/Palma/Downloads/MiddEval3/trainingQ/Teddy/im0.png'); 
%imgL = imread('/Users/Palma/Documents/Valencia/Results/Telescope_v2/EI_3.bmp');
%imgL = imread('/Users/Palma/Documents/Valencia/Results/Box/EI_3.bmp');
%imgC = imread('/Users/Palma/Downloads/MiddEval3/trainingQ/Teddy/im1.png');
%imgC = imread('/Users/Palma/Documents/Valencia/Results/Telescope_v2/EI_4.bmp');
imgR = imread('/Users/Palma/Downloads/MiddEval3/trainingQ/Teddy/im1.png');
%imgR = imread('/Users/Palma/Documents/Valencia/Results/Box/EI_5.bmp');
%imgR = imread('/Users/Palma/Documents/Valencia/Results/Telescope_v2/EI_5.bmp');

imgL01 = double(imgL) / 255;
imgR01 = double(imgR) / 255;
%imgR01 = double(imgC) / 255;
imgL = rgb2gray(imgL);
imgR = rgb2gray(imgR);
%imgC = rgb2gray(imgC);
imgLdenoised = medfilt2(imgL);
imgRdenoised = medfilt2(imgR);

cd DEPTH/CPP/
%mex stereo_matching.cpp
%mex tri_stereo_v2.cpp
%mex bin_stereo.cpp
mex bin_stereo_v2.cpp
mex calc_sad_vol.cpp
%mex copyimg.cpp
cd ../..


%img1 = copyimg(double(imgL), double(imgR));
%imagesc(img1);
cv_v2 = bin_stereo_v2(double(imgL), double(imgR), 0, 60);

%sad_vol = calc_sad_vol(double(imgL), double(imgR), 9);
[~, dcv] = min(cv_v2, [], 3);
dcvgf = zeros(size(cv_v2,1), size(cv_v2,2), size(cv_v2,3));
cv_test = zeros(size(cv_v2,1), size(cv_v2,2), size(cv_v2,3));
for i = 1:size(cv_v2,3)
    imgTranslate = imtranslate(imgR, [(i-1) 0],'FillValues',255);
    a = compute_diff_image_aw(double(imgL), double(imgTranslate), 13, 0.5);
    cv_test(:,:,i) = a;
    dcvgf(:,:,i) = imguidedfilter(cv_v2(:,:,i));
end
[~, dtest] = min(cv_test, [], 3);
[~, dcv2] = min(dcvgf, [], 3);
imagesc(dtest), figure,imagesc(dcv, [0 15]);
sadimg = calc_sad(double(imgL), double(imgR), 9);
nccimg = calc_ncc(double(imgL), double(imgR), 9);
figure, imagesc(sadimg), figure, imagesc(nccimg);

[sad, census] = bin_stereo(double(imgL), double(imgR), 15, 50);
[~, dS] = min(sad, [], 3);
[~, dC] = min(census, [], 3);
figure, imagesc(dC)
figure, imagesc(dS)
[L, R] = tri_stereo_v2(double(imgL), double(imgC), double(imgR), 10, 25, 15, 1, 1);
[~, dL] = min(L, [], 3);
figure, imagesc(dL);
[~, dR] = min(R, [], 3);
figure, imagesc(dR);
[L1, R1] = tri_stereo(double(imgL), double(imgC), double(imgR), 0, 64, 5, 1, 1);
[L2, R2] = tri_stereo(double(imgL), double(imgC), double(imgR), 0, 64, 5, 0, 1);
ddd = disparity(imgL, imgR);
[~, dL] = min(L, [], 3);
[sad, census, ncc, ssd] = stereo_matching(double(imgL), double(imgC), 5, 25);
fprintf('first stereo');
[sad2, census2, ncc2, ssd2] = stereo_matching(double(imgC), double(imgL), -25, -5);
fprintf('second stereo');
[sad3, census3, ncc3, ssd3] = stereo_matching(double(imgC), double(imgR), 5, 25);
fprintf('third stereo');

[c, dispL] = min(0.7.*sad+0.3.*census, [], 3);
[c, disp2] = min(0.7.*sad2+0.3.*census2, [], 3);
dispR = size(sad,3) - disp2;
[c, dispRR] = min(0.7.*sad3+0.3.*census3, [], 3);

figure, imagesc(dispL);
figure, imagesc(dispR);
figure, imagesc(dispRR);
similar_values = (3-(abs(dispR-dispRR))>0 .* (dispR>1));
figure, imagesc(dispR.*similar_values);



[c, disp4] = min(sad, [], 3);
[c, disp5] = min(ssd, [], 3);

figure, imagesc(disp);
figure, imagesc(disp2);
figure, imagesc(disp3);

[sad, census, ncc, ssd] = stereo_matching(double(imgL), double(imgR));

[c, disp] = min(sad, [], 3);
[c, disp2] = min(ssd, [], 3);
[c, disp3] = min(0.7.*sad+0.3.*census, [], 3);
figure, imagesc(disp);
figure, imagesc(disp2);
figure, imagesc(disp3);

%[sad2, census2] = stereo_matching(double(imgLdenoised), double(imgRdenoised));


