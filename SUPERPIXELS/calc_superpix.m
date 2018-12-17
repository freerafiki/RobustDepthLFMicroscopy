function superpixels1 = calc_superpix(im, mask, sp_size, sp_reg, sp_type)
fprintf('\nCalculating superpixels..\n');
tic
% IM contains the image in RGB format as before
%im = im2single(img);
regionSize = sp_size; %sp_size;
regularizer = sp_reg;
%superpixels1 = vl_slic(single(uint8(im)), regionSize, regularizer) ;
%regularizer1 = 0.01;
%regularizer2 = 1.0;
%gfimg = imguidedfilter(uint8(im), 'NeighborhoodSize', [7 7]);
%superpixels1 = vl_slic(single(uint8(gfimg)), regionSize, regularizer1) ;
%sp_num = size(im, 1) * size(im, 2) / (500);

if sp_type == 1

    uint8im = uint8(im);
    gfimg = imguidedfilter(uint8im);
    superpixels1 = vl_slic(single(uint8(gfimg)), regionSize, regularizer) ;
    %[superpixels1, N] = slicmex(gfimg, 500, 20);
    
else
    
    uint8im = uint8(im);
    gfimg = imguidedfilter(uint8im);
    %superpixels1 = vl_slic(single(uint8(gfimg)), regionSize, regularizer) ;
    number_of_sp = 500; %%size(im,1) * size(im,2) / sp_size; % use 500
    compactness_factor = 20;
    [superpixels1, N] = slicmex(gfimg, number_of_sp, compactness_factor);
end


%[L2, N2] = slicomex(uint8im, 500);
%{
subplot(221)
imagesc(uint8(im))
subplot(222)
imagesc(superpixels1)
subplot(223)
imagesc(L)
subplot(224)
imagesc(L)
fprintf('a');
%}
%superpixels1 = superpixels1.*uint32(im2uint8(mask)./255) ;
%regionSize2 = 50;
%regularizer2 = 0.1;
%superpixels2 = vl_slic(single(im), regionSize2, regularizer2) ;
%superpixels2 = superpixels2.*uint32(im2uint8(mask)./255) ;

%subplot(121)
%imagesc(superpixels2)
%subplot(122)
%imagesc(superpixels2.*uint32(im2uint8(mask)./255))
%show the image

%{
imlab = vl_xyz2lab(vl_rgb2xyz(im)) ;
segments = vl_slic(single(imlab), 10, 0.1) ;
%show the iamge

perim = true(size(im,1), size(im,2));
for k = 1 : max(segments(:))
    regionK = segments == k;
    perimK = bwperim(regionK, 8);
    perim(perimK) = false;
end

perim = uint8(cat(3,perim,perim,perim));

finalImage = im .* double(perim);
figure()
imagesc(finalImage);

%}
fprintf('Done!                       ');
toc