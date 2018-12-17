function save_fs(folder_path, fs, map)

for i = 1:size(fs,4)
    name = strcat(folder_path, 'fs_', num2str(i), '_f=', num2str(map(i)), '.bmp');   
    imwrite(uint8(fs(:,:,:,i)),name);
end
