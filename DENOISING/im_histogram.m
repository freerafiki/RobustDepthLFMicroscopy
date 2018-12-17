function [histogram] = im_histogram(input_img)
% for each pixel in the input image
% calculate the histogram for this image

h = size(input_img,1); % height
w = size(input_img,2); % width
histogram = zeros(256, 1);
for r = 1:h
    for c = 1:w
        histogram(input_img(r, c) + 1) = ...
            histogram(input_img(r, c) + 1) + 1;
    end
end

end