function [depth_D, depth_C] = initial_depth_MGM(cvD, cvC, img, matte)

%% TODO

%% PARAMETERS
% MGM_wrapper - calls MGM approximate optimization
%
%   labeling = MGM_wrapper(unary, NDIR, P1, P2, MGM, VTYPE, w)
%
%   MGM solves
%         E(D) = \sum_p C(p,D_p) + \sum_{(p,q)} w(p,q) V (D_p, D_q)
%   which is defined on the 4- or 8-connected image grid of size W x H.
%   D: label map (image), D_p denotes the label of pixel p.
%      The labels are numbered from 0 to L-1
%   C: unary term in the form of a cost volume of size W x H x L. 
%      For each pixel p and each label gives a matching cost.
%   V: the interaction potential with parameters P1 and P2
%      One option is to use SMG's potential (Hirschmuller'08)
%                     | 0  if |a-b|==0
%           Vh(a,b) = |P1  if |a-b|==1
%                     |P2  otherwise
%
%      another option is to use a truncated linear potential (Felzenszwalb-Huttenlocher'06)
%           Vl(a,b) =  min(P1*|a-b|, P2)
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
%   INPUTS: 
%   unary  the cost volume C of size W x H x L
%   NDIR   pass directions 4:
%                          8: (default) 
%   P1,P2  interaction potential params: P1=8, P2=32 (default)
%   MGM    # of messages   1: 1D propagation in each pass as SGM
%                          2: use messages from 2 neighbors (default)
%   VTYPE  V potential     0: Vh as in SGM (Hirschmuller'08) (default)
%                          1: Vl truncated lienar (Felzenszwalb-Huttenlocher'06)
%   w      edge weights of size W x H x 8 (default: all ones)
%
%   OUTPUTS:
%   labeling: approximated MAP solution 
%
%
%   For more details see:  http://dev.ipol.im/~facciolo/mgm
%
%   Copyright (C) 2015, Gabriele Facciolo

NDIR = 8; % Number of directions
P1 = 0.5; % Penalty 1
P2 = 1.0; % Penalty 2
MGM = 4; % 1, 2 or 4 messages
VTYPE = 2; % TLP

w = create_edge_costs(img);

cd /data1/palmieri/MGM
depth_D = MGM_wrapper(cvD, NDIR, P1, P2, MGM, VTYPE, w);
depth_C = MGM_wrapper(cvC, NDIR, P1, P2, MGM, VTYPE, w);
cd /data1/palmieri/Valencia/Code/matlab_depth

%% apply matte?
depth_D = depth_D .* double(matte);
depth_C = depth_C .* double(matte);
