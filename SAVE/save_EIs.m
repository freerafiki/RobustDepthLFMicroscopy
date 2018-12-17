function save_EIs(folder_path, EIs)

for i = 1:size(EIs,4)
    name = strcat(folder_path, 'EI_', num2str(i), '.bmp');   
    imwrite(uint8(EIs(:,:,:,i)),name);
end