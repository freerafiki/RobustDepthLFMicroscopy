function smoothed = smooth_cost_volume(cv, superpixels, img)

    %% SUPERPIXELS/GRID SMOOTHING PART
    tic
    sp_max_val = max(max(superpixels));
    [garbage, disp] = min(cv,[],3);
    gimg = rgb2gray(uint8(img));
    num_of_imgs = size(cv, 3);
    smoothed = zeros(size(cv,1), size(cv,2), size(cv,3));
    for i = 1: sp_max_val
        % get pixels
        indices = superpixels == i;
        img = disp(superpixels==i);
        subplot(221)
        imagesc(double(indices).*double(gimg));
        bins = linspace(0, num_of_imgs, (num_of_imgs/2));
        subplot(222)
        histogram(img, bins);
        [counts, edges] = histcounts(img, bins);
        sorted = sort(counts);
        z = ones(1,num_of_imgs+1);
        if size(counts,2) > 0
            peak1 = sorted(size(counts,2));
            if size(counts,2) > 1
                peak2 = sorted(size(counts,2)-1);
            end
            if size(counts,2) > 2
                peak3 = sorted(size(counts,2)-2);
            end
            ind1 = 0; ind2 = 0; ind3 = 0;
            for j = 1:size(counts, 2)
                if counts(j) == peak1
                    ind1 = j;
                elseif counts(j) == peak2
                    ind2 = j;
                elseif counts(j) == peak3
                    ind3 = j;
                end
            end
            sigma = 3;
            x = linspace(0, num_of_imgs, num_of_imgs+1);
            
            z = z + max(0, 1-((ind1-x)./sigma).^2);
            if peak2 > peak1*.5 && size(counts,2)>1
                z = z + max(0, 1-((ind2-x)./sigma).^2);
            end
            if peak3 > peak1*.5 && size(counts,2)>2
                z = z + max(0, 1-((ind3-x)./sigma).^2);
            end
        end
        %smooth = z;
        subplot(2,2,[3 4])
        plot(z);
        for k = 1:size(cv,3)
            smoothed(:,:,k) = smoothed(:,:,k) + indices.*z(k);
        end
    end

fprintf('Done!          ')
toc
end
