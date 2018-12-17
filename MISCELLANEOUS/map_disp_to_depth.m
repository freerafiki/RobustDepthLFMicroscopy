function dfC_depth = map_disp_to_depth(dfC_disp)

%% Disparity values are projected to depth using the formula 
%% depth = offset - disp*step_pix

step_pix = 13.5;
offset = 70;
mex mapd2d.cpp
dfC_depth = mapd2d(dfC_disp, step_pix, offset);