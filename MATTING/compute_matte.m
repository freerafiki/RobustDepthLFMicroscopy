function matte = compute_matte(aif, trimap)

fprintf('\nCompute the matte..\n');
tic
matte_already_ready = false;
if matte_already_ready == false
    matte = fast_matting(aif, trimap);
    thresh = graythresh(matte);
    matte = im2bw(matte, thresh);
else
    matte = imread('/data1/palmieri/Valencia/Images/Type1/FocalStack3/matte/matte.png');
    matte = im2bw(matte, 0.01);
end
%in case of having the image processing toolbox
%matte = imbinarize(matte, 1); 
%matte = im2uint8(matte);

%uncomment to show trimap
%{

RGB = label2rgb(matte); 	 
figure;
imshow(RGB)
axis off
title('RGB Segmented Image')

%}

fprintf('Done!                       ');
toc