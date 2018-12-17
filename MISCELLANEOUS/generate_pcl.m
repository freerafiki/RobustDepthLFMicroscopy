function generate_pcl(depthmap, image, matte)

for i=1:3
    masked_img(:,:,i) = image(:,:,i).*double(matte);
end
masked_depth = depthmap.*double(matte);
num_of_points = sum(sum(matte>0));
fileID = fopen('/Users/Palma/Documents/Valencia/TestImages/3/xyz.txt','w');
formatSpec = '%3.3f %3.3f %3.3f\n';
[x,y] = meshgrid(1:size(image,1), 1:size(image,2));
%fprintf(fileID,'ply\nformat ascii 1.0\nelement vertex 466489\nproperty float32 x\nproperty float32 y\nproperty float32 z\n');
%fprintf(fileID,'property uchar red\nproperty uchar green\nproperty uchar blue\nelement face 0\n');
%fprintf(fileID,'property list uint8 int32 vertex_index\nend header\n');
fprintf(fileID, formatSpec,x,y,masked_depth);
%fprintf(fileID,'\n');