function matte = fast_matting(aif, trimap)

if size(size(aif),2) == 3
    aif = rgb2gray(uint8(aif));
end
level = graythresh(aif);
matte = im2bw(aif, level);
matte = double(matte);