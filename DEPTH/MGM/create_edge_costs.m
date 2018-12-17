function w = create_edge_costs(img)

% MGM_wrapper - calls MGM approximate optimization
%
%   w: are the graph's edge "weights". For each pair of negihboring 
%      pixels (p,q), w(p,q) gives the corresponding weight.
%      w can actually be used to adapt the parameters P1 and 
%      P2 on the pixel basis like:
%                 V (D_p, D_q, P1(w(p)), P2(w(p)) )
%      But currently it just multiplies the potential.
%      The weights are represented with a stack of 8 images. 
%      For a pixel p each image of the stack correspond the weights to 
%      int neighbouring pixels: W, E, S, N, (NW, NE, SE, SW)
%      For a 4-connectivity, the first 4 images are the weights to the
%      edges, while for 8-connectivity the stack the 8 images are used
%
%                         NEIGHBORS
%
%            (-1,-1) NW    (0,-1)      NE(1,-1) 
%                            |N
%                            |
%            (-1,0) W ---    o    ---  E (1,0)
%                            | 
%                            |S
%            (-1,1) SW     (0,1)       SE(1,1) 
%
%   Copyright (C) 2015, Gabriele Facciolo

width = size(img, 2);
height = size(img, 1);
w = zeros(height, width, 8);
w(:,2:width,1)   = mean( imabsdiff(img(:,2:width,:), img(:,1:width-1,:)),3);
w(:,1:width-1,2) = mean( imabsdiff(img(:,1:width-1,:), img(:,2:width,:)),3);
w(1:height-1,:,3)= mean( imabsdiff(img(1:height-1,:,:), img(2:height,:,:)),3);
w(2:height,:,4)  = mean( imabsdiff(img(2:height,:,:), img(1:height-1,:,:)),3);
w(2:height,2:width,5)  = mean( imabsdiff(img(2:height,2:width,:), img(1:height-1,1:width-1,:)),3);
w(2:height,1:width-1,6)  = mean( imabsdiff(img(2:height,1:width-1,:), img(1:height-1,2:width,:)),3);
w(1:height-1,1:width-1,7)  = mean( imabsdiff(img(1:height-1,1:width-1,:), img(2:height,2:width,:)),3);
w(1:height-1,2:width,8)  = mean( imabsdiff(img(1:height-1,2:width,:), img(2:height,1:width-1,:)),3);
% w...,5)...
% w = exp(-w.*w/64)*4;
