function save_img_and_depth_jpg(folder_path, depth, img, name, matte)

name_depth = strcat(folder_path, name, '_depth.png'); 
name_depth2 = strcat(folder_path, name, '_no_norm_depth.png'); 
name_depth3 = strcat(folder_path, name, '_col_depth.png');
name_depth4 = strcat(folder_path, name, '_col_depth_inv.png'); 
name_img = strcat(folder_path, name, '.png'); 
maxval = max(max(depth));
imgtosave = depth./maxval.*4;

%matted_depth = depth - matte;
%[hist, counts] = histcounts(matted_depth);

imwrite(imgtosave,name_depth);
imwrite(uint8(depth), name_depth2);
imwrite(uint8(img),name_img);
imwrite(uint8(depth.*2), colormap(jet), name_depth3);

inv_depth = 50 - depth;
imwrite(uint8(inv_depth), colormap(jet), name_depth4);