function save_depth_CV(folder_path, cost_volume, name, matte)

for i = 1:size(cost_volume,3)
    cost_volume(:,:,i) = imguidedfilter(cost_volume(:,:,i));
end
[~, depth] = min(cost_volume, [], 3);
name1 = strcat(folder_path, name, '.png');  
name2 = strcat(folder_path, name, '_no_norm.png');  
%name3 = strcat(folder_path, name, '_MGM.png');  

%{
%% PARAMETERS
NDIR = 8; % Number of directions
P1 = 0.5; % Penalty 1
P2 = 1.0; % Penalty 2
MGM = 4; % 1, 2 or 4 messages
VTYPE = 2; % TLP
w = create_edge_costs(depth);
cd /data1/palmieri/MGM
depthmap = MGM_wrapper(cost_volume, NDIR, P1, P2, MGM, VTYPE, w);
cd /data1/palmieri/Valencia/Code/matlab_depth
%}

maxval = max(max(depth));
imgtosave = depth.*double(matte)./maxval;
imwrite(uint8(imgtosave), colormap(jet), name1);
imwrite(uint8(depth), colormap(jet) ,name2);
%imwrite(depthmap, colormap(jet) ,name3);