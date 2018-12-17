addpath(genpath('DENOISING'));

noisy_img = rgb2gray(imread('/Users/Palma/Documents/Valencia/TestImages/facce/EI_4.bmp'));

medfilt = medfilt2(noisy_img);
Kaverage = filter2(fspecial('average',3),noisy_img)/255;
% non local means
%sigma = 10;
%nl = NLmeansfilter(double(noisy_img),5,2,sigma);
% Set bilateral filter parameters.
w     = 5;       % bilateral filter half-width
sigma = [3 0.1]; % bilateral filter standard deviations
nnn = double(noisy_img) / 255;
bilat = bfilter2(nnn,w,sigma);

imagesc(noisy_img);
figure, imagesc(medfilt);
figure, imagesc(Kaverage);
figure, imagesc(bilat);