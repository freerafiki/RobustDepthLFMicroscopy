function save_matte(folder_path, matte)

name = strcat(folder_path, 'matte.bmp');  
% remap [0,1] to [0,255] for visualization
imwrite(uint8(matte*255),name);