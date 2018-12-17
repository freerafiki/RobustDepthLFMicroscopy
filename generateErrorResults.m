%% Calculate errors of the different methods against ground truth
clear all; close all; clc

%% STEP 1: prepare strings for images, read matte

folder = '/data1/palmieri/2018/November/NumAnalysis/Images';
nameGT = '_gt.png';
namePROP_sm = '_sm.png';
namePROP_df = '_df.png';
namePROP = '_final.png';
nameCNN = '_cnn_dmax29.bin';
nameSFF = '_sff.png';
nameMatte = '_matte.bmp';
numOfImages = 12;
h = 683;
w = 683;

%% 1b) matte

%matte = imtranslate(matte, [-1, -1]);

%% 1c) ERROR IMAGES
err_img_prop_sm = zeros(h, w, numOfImages);
err_img_prop_df = zeros(h, w, numOfImages);
err_img_prop = zeros(h, w, numOfImages);
err_img_cnn = zeros(h, w, numOfImages);
err_img_sff = zeros(h, w, numOfImages);

%% 1d) CALCULATE ERRORS
err_prop_sm = zeros(40704, numOfImages); % 9 - 79953 % 8 - 41665 % 10 - 97113 % 11 - 58605 % 12 - 40704
err_prop_df = zeros(40704, numOfImages);
err_prop = zeros(40704, numOfImages);
err_cnn = zeros(40704, numOfImages);
err_sff = zeros(40704, numOfImages);

%% STEP 2: read image and calculate errors
for j = 12:numOfImages
    
    fprintf('Analyzing image %d/%d..\n', j, numOfImages); 
    % GT HAS TO BE UNDISTORTED
    gtpath = (strcat(folder, filesep, num2str(j), nameGT));
    debugging = false;
    disp = true;
    gt = read_and_undistort_gt(gtpath, debugging, disp);
    gt(gt > 100) = 100;
    gt = double(gt);
    % Proposed method have index inverted (depth/disparity)
    prop_sm = imread(strcat(folder, filesep, num2str(j), namePROP_sm));
    prop_sm = double(prop_sm);
    %prop_sm = double(max(max(prop_sm)) - prop_sm);
    prop_df = imread(strcat(folder, filesep, num2str(j), namePROP_df));
    prop_df = double(prop_df);
    %prop_df = double(max(max(prop_df)) - prop_df);
    prop = imread(strcat(folder, filesep, num2str(j), namePROP));
    prop = double(prop);
    %prop = double(max(max(prop)) - prop);
    % CNN is saved into a binary file that need to be reshaped
    cnnStruct = memmapfile(strcat(folder, filesep, num2str(j), nameCNN), 'Format', 'single');
    cnn = reshape(cnnStruct.Data, [h w])';
    cnn = double(cnn);
    % SFF as proposed needs to invert the values
    sff = imread(strcat(folder, filesep, num2str(j), nameSFF));
    sff = double(max(max(sff)) - sff);
    % read the mask
    matte = imread(strcat(folder, filesep, num2str(j), nameMatte)) ./ 255;
    matte = imerode(matte, [0 1 0; 1 1 1; 0 1 0]);
    matte = double(matte);
    
    %% Fill the error maps and arrays
    err_img_prop_sm(:,:,j) = imabsdiff(prop_sm .* matte, gt .* matte);
    err_img_prop_df(:,:,j) = imabsdiff(prop_df .* matte, gt .* matte);
    err_img_prop(:,:,j) = imabsdiff(prop .* matte, gt .* matte);
    err_img_cnn(:,:,j) = imabsdiff(cnn .* matte, gt .* matte);
    err_img_sff(:,:,j) = imabsdiff(sff .* matte, gt .* matte);
    
    [row, col] = find(matte);
    arraySize = size(row,1);
    err_prop_sm(:,j) = zeros([arraySize 1]);
    err_prop_df(:,j) = zeros([arraySize 1]);
    err_prop(:,j) = zeros([arraySize 1]);
    err_cnn(:,j) = zeros([arraySize 1]);
    err_sff(:,j) = zeros([arraySize 1]);
    
    for l = 1:size(row,1)
        
        err_prop_sm(l,j) = err_img_prop_sm(row(l), col(l), j);
        err_prop_df(l,j) = err_img_prop_df(row(l), col(l), j);
        err_prop(l,j) = err_img_prop(row(l), col(l), j);
        err_cnn(l,j) = err_img_cnn(row(l), col(l), j);
        err_sff(l,j) = err_img_sff(row(l), col(l), j);
    end
    
    nameFile = strcat('/data1/palmieri/2018/November/NumAnalysis/Images/', num2str(j), '_errors.txt');
    fileID = fopen(nameFile,'w');
    fprintf(fileID, 'SM: %f + %f\n', mean(err_prop_sm(:,j)), std(err_prop_sm(:,j)));
    fprintf(fileID, 'DF: %f + %f\n', mean(err_prop_df(:,j)), std(err_prop_df(:,j)));
    fprintf(fileID, 'PROP: %f + %f\n', mean(err_prop(:,j)), std(err_prop(:,j)));
    fprintf(fileID, 'CNN: %f + %f\n', mean(err_cnn(:,j)), std(err_cnn(:,j)));
    fprintf(fileID, 'SFF: %f + %f\n', mean(err_sff(:,j)), std(err_sff(:,j)));
    
    %colormap jet
    subplot(231); imagesc(gt .* matte, [0 19]);  
    subplot(232); imagesc(prop .* matte, [0 19]); title('Prop'); 
    subplot(233); imagesc(cnn .* matte, [0 19]); title('cnn'); 
    subplot(234); imagesc(sff .* matte, [0 19]); title('sff'); 
    subplot(235); imagesc(prop_df .* matte, [0 19]); title('def'); 
    subplot(236); imagesc(prop_sm .* matte, [0 19]); title('sm'); 
    
end

    

