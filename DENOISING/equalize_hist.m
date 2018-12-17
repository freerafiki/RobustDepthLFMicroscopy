function [output_img] = equalize_hist(input_img)
    % The function applies histogram equalization on a gray scale image
    % returning a gray scale image whose histogram is approximately flat.

    L = 256; % number of grey levels used
    h = size(input_img,1); % height of input&output
    w = size(input_img,2); % width of input&output
    number_of_pixels = h * w;

    % create output image that is the same size and data type as input
    output_img = zeros(h, w);
    output_img = cast(output_img, 'like', input_img);

    % the index from 1 to 256 corresponds to intensity from 0 to 255
    cdf = zeros(256, 1);
    equalized = zeros(256, 1);
    histogram = im_histogram(input_img);

    % for each pixel in the input image
    % calculate the cdf for this image
    count = 0;
    for i = 1:size(histogram)
        count = count + histogram(i);
        cdf(i) = count / number_of_pixels;
        equalized(i) = round(cdf(i) * (L - 1));
    end

    % for each pixel in the output image
    % calculate the histogram equalization result
    for r = 1:h
        for c = 1:w
            output_img(r, c) = equalized(input_img(r, c) + 1);
        end
    end


end