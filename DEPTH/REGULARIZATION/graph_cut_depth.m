%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018.06
% Luca Palmieri
% 
% Graph Cuts optimization using gco implementation from
% http://vision.csd.uwo.ca/code/
%
% Smoothness cost taken from 
% 2015.05.12 Hae-Gon Jeon
% Accurate Depth Map Estimation from a Lenslet Light Field Camera
% CVPR 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function depthmap = graph_cut_depth(cost_volume, c_img, param)

    [height, width, num_labels] = size(cost_volume);
    Data = reshape(cost_volume,[],num_labels)';
    h = GCO_Create(height*width, num_labels);  
    
    %% DATA COST 
    % data cost is in int32 so no decimal. We multiply per 4 to obtain a
    % integer number every 0.25 interval
    %dec_interval = 0.25;
    %factor = 1/dec_interval;
    %int32_Data = round(Data .* factor);
    %GCO_SetDataCost(h, int32(int32_Data));
    
    numQuantiz =5*num_labels;
    QQ = ( Quantiz(Data(:) , linspace(min(Data(:)),max(Data(:)),numQuantiz)) ) ;
    Data_idx = reshape(QQ,size(Data));
    GCO_SetDataCost(h, int32(Data_idx .* param.data));
    
    %% SMOOTHNESS COST
    % just every pixel it decreases by one
    max_cost = 10;
    Smoothness = min(triu(gallery('circul', [0:num_labels-1]'),0)+triu(gallery('circul', [0:num_labels-1]'),0)',max_cost);
    GCO_SetSmoothCost(h,int32(Smoothness .* param.smooth));
    
    %% NEIGHBOURS
    % a sparse matrix that contains 1 for neighbour
    % [ 0 1 0 1 0 ]
    % [ 0 0 1 0 1 ]
    % [ 0 0 0 1 0 ]
    % [ 0 0 0 0 1 ]
    % [ 0 0 0 0 0 ]
    % Warning: Neighbours array should be upper-triangular; entries below the diagnonal will be ignored 
    size_hw = height * width;
    tmp_1 = linspace(1, size_hw-1, size_hw-1);
    tmp_2 = linspace(1, size_hw-1, size_hw-1)+1;
    tmp_3 = linspace(1, size_hw-width, size_hw-width);
    tmp_4 = linspace(1, size_hw-width, size_hw-width)+width;
    i = [tmp_1, tmp_3];
    j = [tmp_2, tmp_4];
    v = ones(size(i));
    Neighbours = sparse(i,j,v,size_hw, size_hw);
    GCO_SetNeighbors(h,Neighbours);
    
    %% Run expansion algorithm to get optimal label
    ExpansionResult = GCO_Expansion(h);
    Labels = GCO_GetLabeling(h); %% Labels is the depth map
    [Energy, DataTerm, SmoothnessTerm] = GCO_ComputeEnergy(h);
    GCO_Delete(h);
    
    depthmap = reshape(Labels,[height width]);
    
    
    
    
    
    
    
    