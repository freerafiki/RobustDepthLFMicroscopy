function save_trimap(folder_path, trimap)

name = strcat(folder_path, 'trimap.bmp');  
% remap 1,2,3 to 85, 170, 255 for visualization
imwrite(uint8(trimap*85),name);