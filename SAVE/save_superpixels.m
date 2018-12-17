function save_superpixels(folder_path, superpixels, im)


name = strcat(folder_path, 'superpixels_colored.bmp');  
name2 = strcat(folder_path, 'image_superpixels.bmp');  
perim = true(size(im,1), size(im,2));
for k = 1 : max(superpixels(:))
    regionK = superpixels == k;
    perimK = bwperim(regionK, 8);
    perim(perimK) = false;
end

perim = uint8(cat(3,perim,perim,perim));
finalImage = im .* double(perim);
%imagesc(finalImage);
%imwrite(uin(superpixels),name);
imwrite((finalImage),name2);
imwrite(uint8(superpixels), name);