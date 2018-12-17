%%% Image Denoising using Bilateral Filter and Detail Wavelet Thresholding.
%%% Author : B.K. SHREYAMSHA KUMAR 

%%% Copyright (c) 2012 B. K. Shreyamsha Kumar 
%%% All rights reserved.
 
%%% Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, 
%%% modify, and distribute this code (the source files) and its documentation for any purpose, provided that the 
%%% copyright notice in its entirety appear in all copies of this code, and the original source of this code, 
%%% This should be acknowledged in any publication that reports research using this code. The research is to be 
%%% cited in the bibliography as:

%%% B. K. Shreyamsha Kumar, ?Image Denoising based on Gaussian/Bilateral Filter and its Method Noise Thresholding", 
%%% Signal, Image and Video Processing, pp. 1-14, 2012. (doi: 10.1007/s11760-012-0372-7)


clear all;
close all;
clc;

%%% Bilateral Filter Parameters.
sigmas=1.8;  %%% Spatial (Geometric) Sigma.
mult_factor_sigmar=5; %%% Multiplication factor for SigmaR.
ksize=11;   %%% Kernal Size.

%%% Wavelet Transform Parameters.
Nlevels=3; %% 3 levels.
NoOfBands=3*Nlevels+1;
wname='db8'; %% db8 sym8 db16 coif5 bior6.8
sorh='s'; % s or h or t -> trimmed

x = (imread('/Users/Palma/Documents/Valencia/TestImages/facce/EI_4.bmp'));
% x = imread('Barbara256.png');
% x = imread('boat256.jpg');
% x = imread('baboon256.jpg');

if(size(x,3)==3)
    x=rgb2gray(x);
end
[M,N]=size(x);

%%% Gaussian Noise addition.
'gaussian noise';
mean_val=0;
noise_std=10; %% 10,20,30,40,50.
sizeA = size(x);
randn('seed',662023); %%% Results depend on 'seed' of the random noise.
xn = double(x) + (noise_std*randn(sizeA)) + mean_val;
xn = max(0,min(xn,255));

%%% PSNR Computation of Noisy Image.
xn_mse=sum(sum((double(x)-double(xn)).^2))/(M*N);
xn_psnr=10*log10(255^2./xn_mse);

tic
%%% Noise Level Estimation using Robust Median Estimator.
[ca,ch,cv,cd]=dwt2(xn,wname);
tt1=cd(:)';
median_hh2=median(abs(tt1)); %% HH1->Subband containing finest level diagonal details.
std_dev2=(median_hh2/0.6745);

%%% Bilateral Filtering.
sigmar=mult_factor_sigmar*std_dev2;
% tic
[yb yg]= bilateral_filt2D(xn,sigmas,sigmar,ksize);
% toc
yd=double(xn)-yb;  %%% Bilateral Method Noise.

%%% General Wavelet Decomposition.
dwtmode('per');
[C,S]=wavedec2(yd,Nlevels,wname);
k=NoOfBands;
CW{k}=reshape(C(1:S(1,1)*S(1,2)),S(1,1),S(1,2));
k=k-1;
st_pt=S(1,1)*S(1,2);
for i=2:size(S,1)-1
   slen=S(i,1)*S(i,2);
   CW{k}=reshape(C(st_pt+slen+1:st_pt+2*slen),S(i,1),S(i,2));  %% Vertical
   CW{k-1}=reshape(C(st_pt+1:st_pt+slen),S(i,1),S(i,2));   %% Horizontal
   CW{k-2}=reshape(C(st_pt+2*slen+1:st_pt+3*slen),S(i,1),S(i,2));  %% Diagonal
   st_pt=st_pt+3*slen;
   k=k-3;
end

%%%% BayesShrink Technique.
tt2=CW{1}(:)';
median_hh2=median(abs(tt2)); %% HH1->Subband containing finest level diagonal details.
std_dev2=(median_hh2/0.6745);
cw_noise_var=std_dev2^2; %% var=std^2

for i=1:NoOfBands-1
    thr=bayesthf(CW{i},cw_noise_var);
    yw{i}=threshf(CW{i},sorh,thr,2);    
end
yw{i+1}=CW{i+1};

%%% General Wavelet Reconstruction.
if(isequal(wname,'DCHWT'))
   %%% Discrete Cosine Harmonic Wavelet Reconstruction
   ydr=dchwtf2(yw,-Nlevels);
else
   k=NoOfBands;
   xrtemp=reshape(yw{k},1,S(1,1)*S(1,2));
   k=k-1;
   for i=2:size(S,1)-1
       xrtemp=[xrtemp reshape(yw{k-1},1,S(i,1)*S(i,2)) reshape(yw{k},1,S(i,1)*S(i,2)) reshape(yw{k-2},1,S(i,1)*S(i,2))];
       k=k-3;
   end
   ydr=(waverec2(xrtemp,S,wname));
end
ybw=yb+ydr;
ybwn=uint8(ybw);
yb1=yb;
yb=uint8(yb);

toc 

%%%% MSE Computation.
yb_mse=sum(sum((double(x)-double(yb)).^2))/(M*N);
ybw_mse=sum(sum((double(x)-double(ybwn)).^2))/(M*N);

%%%% PSNR Computation.
wname,noise_std,xn_psnr
psnr_bfilt=10*log10(255^2./yb_mse);
psnr_bf_mnt=10*log10(255^2./ybw_mse);

%%%% Image Quality Index (IQI) Computation.
mean_org=mean(x(:));
var_org=sum((x(:)-mean_org).^2)/(M*N-1);

%%%% IQI for Bilateral Filter.
mean_bfilt=mean(yb(:));
var_bfilt=sum((yb(:)-mean_bfilt).^2)/(M*N-1);
cross_var_bfilt=sum((x(:)-mean_org).*(yb(:)-mean_bfilt))/(M*N-1);
IQI_bfilt=(4*cross_var_bfilt*mean_org*mean_bfilt)/((var_org+var_bfilt)*(mean_org^2+mean_bfilt^2));

%%%% IQI for Bilateral Filter & Wavelet Thresholding.
mean_bwt=mean(ybwn(:));
var_bwt=sum((ybwn(:)-mean_bwt).^2)/(M*N-1);
cross_var_bwt=sum((x(:)-mean_org).*(ybwn(:)-mean_bwt))/(M*N-1);
IQI_bf_mnt=(4*cross_var_bwt*mean_org*mean_bwt)/((var_org+var_bwt)*(mean_org^2+mean_bwt^2));

%%%% Visual Fidelity Index (VIF) Computation.
%vif_bfilt = vifvec(double(x),yb1);
%vif_bf_mnt = vifvec(double(x),ybw);
% vif_bf_mnt = vifvec(double(x),ybwn) %%% uint8

% vif_bfilt,psnr_bfilt,IQI_bfilt
%vif_bf_mnt,psnr_bf_mnt,IQI_bf_mnt

%%% Denoised Plots
% figure,imshow(uint8(xn)),colormap gray
figure,imshow(yb),colormap gray,title('Bilateral Filtering');
figure,imshow(ybwn),colormap gray,title('Bilateral Filtering & Detail Thresholding');

% %%%% Method Noise Plots;
mn=xn-double(x);
mnb=yd;
mnbwn=xn-double(ybwn);

% figure,imshow(mnb),colormap gray,title('Method Noise-Bilateral Filtering');

% figure,imshow(mnbwn),colormap gray,title('Method Noise-Bilateral Filtering & Detail Thresholding');

