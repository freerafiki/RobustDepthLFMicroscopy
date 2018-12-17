%% %% MANISH KUMAR SHARMA, GCET Greater Noida, 2014 %% weighted median filter on image with salt & pepper noise %% I = rgb2gray(imread('lena.jpg'));
I = im2double(final_depth_FILL);
Z = I; %imnoise(I,'salt & pepper',0.02); % adding Noise 
mask_size = 3;
mask_cent = ceil(mask_size/2);
mask_start = mask_cent - 1;
W = ones(mask_size);
W(mask_cent, mask_cent) = 4;
sum_W = sum(W(:));
W = W./sum_W; % the values should taken like that, the total sum of values of filter % is divided by there sum and value should be equal to 1.
image_padded = padarray(Z,[mask_start, mask_start],'symmetric','both');
[row, col] = size(Z);
row_pad = row + 2 * mask_start;
col_pad = col + 2 * mask_start;
b = Z;
for x = mask_cent:row_pad - mask_start
for y = mask_cent:col_pad - mask_start
%% To make a 3x3 weighted mask into a 1x9 mask
patch_selected = image_padded(x - mask_start: x + mask_start, y - mask_start: y + mask_start);
a1 = W.*patch_selected;
med = median(a1(:))*sum_W; % the5th value is the weighted median
b(x,y) = med;
end
end
figure(1); imshow((Z))
figure(2); imshow(b) %on image with salt & pepper noise %%
%I = rgb2gray(imread('lena.jpg'));
I = im2double(final_depth_MRF);
Z = imnoise(I,'salt & pepper',0.02); % adding Noise
mask_size = 3;
mask_cent = ceil(mask_size/2);
mask_start = mask_cent - 1;
W = ones(mask_size);
W(mask_cent, mask_cent) = 4;
sum_W = sum(W(:));
W = W./sum_W; % the values should taken like that, the total sum of values of filter % is divided by there sum and value should be equal to 1.
image_padded = padarray(Z,[mask_start, mask_start],'symmetric','both');
[row, col] = size(Z);
row_pad = row + 2 * mask_start;
col_pad = col + 2 * mask_start;
b = Z;
for x = mask_cent:row_pad - mask_start
for y = mask_cent:col_pad - mask_start
%% To make a 3x3 weighted mask into a 1x9 mask
patch_selected = image_padded(x - mask_start: x + mask_start, y - mask_start: y + mask_start);
a1 = W.*patch_selected;
med = median(a1(:))*sum_W; % the5th value is the weighted median b(x,y) = med;
end
end
figure(1); imshow((Z))
figure(2); imshow(b)