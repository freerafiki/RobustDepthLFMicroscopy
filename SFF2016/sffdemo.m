clear, clc, close all
%% Load the variables of the focus sequence:
load imdata.mat

% This loads a structure with two fields: images, a 1x30 cell
% array where each cell is a string with the path of one
% frame of the focus sequence; and focus, a 1x30 vector with
% the focus position (in meters) corresponding to each frame
% of the focus sequence. This sequence was generated using:
% http://www.mathworks.com/matlabcentral/fileexchange/55095-defocus-simulation

%% Preview the image sequence:
showimages(imdata.images, imdata.focus, imdata.ROI);

%% Perform SFF using 3-point gaussian interpolation
% as originally proposed by [1] and compute
% reliability according to [2]:

[z, r] = sff(imdata.images, 'focus', imdata.focus);

%% Carve depthmap by removing pixels with R<20 dB:
zc = z;
zc(r<20) = nan;

%% Display the result:
close(gcf), figure
subplot(1,2,1), surf(z), shading flat, colormap copper
set(gca, 'zdir', 'reverse', 'xtick', [], 'ytick', [])
axis tight, grid off, box on
zlabel('pixel depth (mm)')
title('Full depthmap')

subplot(1,2,2), surf(zc), shading flat, colormap copper
set(gca, 'zdir', 'reverse', 'xtick', [], 'ytick', [])
axis tight, grid off, box on
zlabel('pixel depth (mm)')
title('Carved depthmap (R<20dB)')


% References:
% [1] S.K. Nayar and Y. Nakagawa, PAMI, 16(8):824-831, 1994. 
% [2] S. Pertuz, D. Puig, M. A. Garcia, Reliability measure for shape-from
% focus, IMAVIS, 31:725�734, 2013.