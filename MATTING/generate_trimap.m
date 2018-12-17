function trimap = generate_trimap(aif)

fprintf('\nComputing the trimap..\n')
% create trimap from scene
lab_img = rgb2lab(aif);
L = lab_img(:,:,1);
thresh = multithresh(L,2);
seg_I = imquantize(L,thresh);

%uncomment to show trimap
%{

RGB = label2rgb(seg_I); 	 
figure;
imshow(RGB)
axis off
title('RGB Segmented Image')

%}

trimap = seg_I;
fprintf('Done!                       ');
toc