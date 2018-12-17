function depthmap = compute_depth_with_MGM(cv1, cv2, img, matte, D, WS)

%% USE MGM TO COMPUTE DEPTH MAP

%% PARAMETERS
NDIR = 8; % Number of directions
P1 = 0.5; % Penalty 1
P2 = 1.0; % Penalty 2
MGM = 4; % 1, 2 or 4 messages
VTYPE = 2; % TLP

final_cost_volume = merge_two_cv(cv1, cv2, matte, D, WS);
w = create_edge_costs(img);

cd /data1/palmieri/MGM
depthmap = MGM_wrapper(final_cost_volume, NDIR, P1, P2, MGM, VTYPE, w);
cd /data1/palmieri/Valencia/Code/matlab_depth