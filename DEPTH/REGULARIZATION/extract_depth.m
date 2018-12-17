function final_depth = extract_depth(cost_volume, matte, c_img, method)


%% EXTRACT DEPTH FROM COST VOLUME

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       %
%  method:                              %
%  'wta' = Winner Takes All             %
%  'mgm' = More Global Matching [1]     %
%  'gc' = Graph Cuts [2,3]              %
%                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [1] MGM: A Significantly More Global Matching for Stereovision, Gabriele Facciolo, Carlo de Franchis and Enric Meinhardt
%       code at: https://github.com/gfacciol/mgm
% [2] Code from: http://vision.csd.uwo.ca/code/
% [3] Accurate Depth Map Estimation from a Lenslet Light Field Camera, Jeon et al.
%       code at: https://sites.google.com/site/hgjeoncv/home/depthfromlf_cvpr15

if strcmp(method,'wta')
    
    [~, depth] = min(cost_volume, [], 3);
    final_depth = depth .* double(matte);
    
elseif strcmp(method,'mgm')
    
    NDIR = 8; % Number of directions
    P1 = 0.1; % Penalty 1
    P2 = 0.5; % Penalty 2
    MGM = 4; % 1, 2 or 4 messages
    VTYPE = 2; % TLP - Truncated Linear Potential

    w = create_edge_costs(c_img);
    depth = MGM_wrapper(cost_volume, NDIR, P1, P2, MGM, VTYPE, w);
    final_depth = depth .* double(matte);
    
elseif strcmp(method,'gc')


    param.data = 1.0; %0.4;
    param.smooth = 3; %4;
    param.neigh = 0.015;
    %depth = GraphCuts(cost_volume, c_img, param);
    depth = graph_cut_depth(cost_volume, c_img, param);
    final_depth = depth .* int32(matte);
    
else
    
    fprintf('Unknown method!')
    depth = zeros(size(c_img, 1), size(c_img,2));
    final_depth = depth .* double(matte);
    
end
