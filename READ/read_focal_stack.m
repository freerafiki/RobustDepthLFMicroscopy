function [fs, aif] = read_focal_stack(path)

for i = 1:7
    name = strcat(path,int2str(i),'.jpg');
    fs(:,:,:,i) = imread(name);
end
aifname = strcat(path,'beers.jpg');
aif = imread(aifname);