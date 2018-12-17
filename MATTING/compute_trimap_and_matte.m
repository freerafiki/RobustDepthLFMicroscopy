function [trimap, matte] = compute_trimap_and_matte(img, img_type, pitch, matting)

if (img_type == 1)
    % we need the matting
    if (matting == 1)
        % use the three-layer matting
        trimap = generate_trimap(img);
        matte = three_layer_matting(img, trimap);
        thresh = graythresh(matte);
        matte = im2bw(matte, thresh);
    else
        % fast matting
        trimap = generate_trimap(img);
        matte = compute_matte(img, trimap);
    end
elseif (img_type == 3)
    % we don't need the matte, just the circular
    [X,Y] = meshgrid(1:size(img,1), 1:size(img,2));
    X2 = X - size(img,1) / 2;
    Y2 = Y - size(img,2) / 2;
    radius = floor(473 / 2); % for type 3!
    matte = (sqrt(X2.^2+Y2.^2)<radius);
    trimap = matte;
elseif (img_type == 4)
    %trimap = generate_trimap(img);
    %matte = compute_matte(img, trimap);
    [X,Y] = meshgrid(1:size(img,1), 1:size(img,2));
    X2 = X - size(img,1) / 2;
    Y2 = Y - size(img,2) / 2;
    radius = floor(687 / 2); % for type 4!
    matte = (sqrt(X2.^2+Y2.^2)<radius);
    trimap = matte;
elseif (img_type == 5)
    % we don't need the matte, just the circular
    [X,Y] = meshgrid(1:size(img,1), 1:size(img,2));
    X2 = X - size(img,1) / 2;
    Y2 = Y - size(img,2) / 2;
    radius = floor(420 / 2); % for type 5!
    matte = (sqrt(X2.^2+Y2.^2)<radius);
    trimap = matte;
elseif (img_type == 6)
    % we don't need the matte, just the circular
    [X,Y] = meshgrid(1:size(img,1), 1:size(img,2));
    X2 = X - size(img,1) / 2;
    Y2 = Y - size(img,2) / 2;
    radius = floor(pitch / 2); % for type 5!
    matte = (sqrt(X2.^2+Y2.^2)<radius);
    trimap = matte;
elseif (img_type == 7)
    % we don't need the matte, just the circular
    [X,Y] = meshgrid(1:size(img,1), 1:size(img,2));
    X2 = X - size(img,1) / 2;
    Y2 = Y - size(img,2) / 2;
    radius = floor(pitch / 2); % for type 5!
    matte = (sqrt(X2.^2+Y2.^2)<radius);
    trimap = matte;
elseif (img_type == 8)
    % we don't need the matte, just the circular
    [X,Y] = meshgrid(1:size(img,1), 1:size(img,2));
    X2 = X - size(img,1) / 2;
    Y2 = Y - size(img,2) / 2;
    radius = floor(pitch / 2);
    matte = (sqrt(X2.^2+Y2.^2)<radius-10);
    trimap = matte;
elseif (img_type == 9)
    % we need the map
    if (matting == 1)
        % use the three-layer matting
        trimap = generate_trimap(img);
        matte = three_layer_matting(img, trimap);
        thresh = graythresh(matte);
        matte = im2bw(matte, thresh);
    else
        % fast matting
        trimap = generate_trimap(img);
        matte = compute_matte(img, trimap);
    end
elseif (img_type == 10)
    % we need the matting
    if (matting == 1)
        % use the three-layer matting
        trimap = generate_trimap(img);
        matte = three_layer_matting(img, trimap);
        thresh = graythresh(matte);
        matte = im2bw(matte, thresh);
    else
        % fast matting
        trimap = generate_trimap(img);
        matte = compute_matte(img, trimap);
    end
elseif (img_type == 11)
    % we need the matting
    if (matting == 1)
        % use the three-layer matting
        trimap = generate_trimap(img);
        matte = three_layer_matting(img, trimap);
        thresh = graythresh(matte);
        matte = im2bw(matte, thresh);
    else
        % fast matting
        trimap = generate_trimap(img);
        matte = compute_matte(img, trimap);
    end
end