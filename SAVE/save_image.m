function save_image(folder_path, fs, map, EIs)

for i = 1:size(fs,4)
    name = strcat(folder_path, 'fs_', num2str(i), '_f=', num2str(map(i)), '.png');   
    imwrite(uint8(fs(:,:,:,i)),name);
end

for i = 1:size(EIs,4)
    name = strcat(folder_path, 'EI_', num2str(i), '.png');   
    imwrite(uint8(EIs(:,:,:,i)),name);
end