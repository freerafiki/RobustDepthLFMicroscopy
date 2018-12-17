clear
clc
clf
colormap(gray)

% create example image
ima=double(imread('/Users/Palma/Documents/Valencia/TestImages/facce/EI_4.bmp'));
ima(50:100,:)=50;
ima(:,50:100)=2*ima(:,50:100);
fs=fspecial('average');
ima=imfilter(ima,fs,'symmetric');

% add some noise
sigma=10;
rima=ima+sigma*randn(size(ima));

% show it
imagesc(rima)
drawnow

% denoise it
fima=NLmeansfilter(ima,5,2,sigma);

% show results
clf
subplot(2,2,1),imagesc(ima),title('original');
subplot(2,2,2),imagesc(rima),title('noisy');
subplot(2,2,3),imagesc(fima),title('filtered');
subplot(2,2,4),imagesc(rima-fima),title('residuals');

