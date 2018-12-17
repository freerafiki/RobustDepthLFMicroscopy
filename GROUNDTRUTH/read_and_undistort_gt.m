function undistorted_gt = read_and_undistort_gt(gt_path, debug, disp)

gt = imread(gt_path);

% convert to rgb
%gtgray = rgb2gray(gt);
%gtgray = gtgray - min(min(gtgray));
gtgray = gt(:,:,1);
[h,w] = size(gtgray);
gtundist = zeros(size(gtgray));
gtdist = zeros(size(gtgray));
amax = 68.666 / 2; % field of view
pmax = floor(w/2); % max pixel distance from center
% angle(pmax) = amax
center = ceil(w/2);

% correct distortion
for i=1:w
    for j=1:h
        
        pdist = sqrt( power(i-center, 2) + power(j-center,2) ); % pixels distance
        if pdist > pmax
            gtundist(i,j) = 0.0;
            gtdist(i,j) = 0.0;
        else
            
            %angle
            angle = pdist / pmax * amax;

            % correct value
            gtundist(i,j) = gtgray(i,j) * cos(degtorad(angle));
            gtdist(i,j) = gtgray(i,j);
        end
        %if (mod(i,100) == 0) && (mod(j,100) == 0)
        %   sprintf('index %d, %d, angle %.3f, cos %.3f, gt %.3f, corrected %.3f', i,j, angle, cos(degtorad(angle)), gtgray(i,j), gtundist(i,j))
        %end
    end
end

if debug
    figure, imagesc(gtgray), figure, imagesc(gtundist)
end

% convert to float
depth = double(gtgray);
% normalize to 0-1
depthnorm = depth / max(max(depth));
% get to 0 -255
depthnorm255 = depthnorm * 255 + 1;
if ~disp
    undistorted_gt = 1900 ./ gtundist ;
else
    undistorted_gt = gtundist ./10;
end
if debug
    figure, mesh(gtdist), colormap jet 
    figure, mesh((gtundist)), colormap jet
end