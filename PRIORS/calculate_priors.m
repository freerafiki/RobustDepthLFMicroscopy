function priors = calculate_priors(refocused_central_img, matte, maptype)

%% CALCULATE PRIORS ON IMAGE

if length(size(refocused_central_img)) == 3
    img = rgb2gray(uint8(refocused_central_img));
else
    img = refocused_central_img;
end

%% OVER AND UNDEREXPOSED PARTS
underexposed_mask = double(img > 10); % usually the matte
overexposed_mask = double(img < 210);

%% EDGE MAP
if maptype == 0
    edgemap = ones(size(img, 1), size(img, 2));
elseif maptype == 1
    % DoG
    w_size = 11;
    sigma_1 = 10;
    sigma_2 = 20;
    H1 = fspecial('gaussian',w_size,sigma_1);
    H2 = fspecial('gaussian',w_size,sigma_2);
    dog = H1 - H2;
    edgemap = conv2(double(img),dog,'same');
    %imagesc(dogFilterImage1)
end

combined_map = edgemap; % .* underexposed_mask .* overexposed_mask;

priors = abs(combined_map) .* double(matte);