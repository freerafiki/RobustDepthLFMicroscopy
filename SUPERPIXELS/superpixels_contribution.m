function cv_sp = superpixels_contribution(cv, superpixels, img, matte, sigma);

    %% SUPERPIXELS/GRID SMOOTHING PART
    fprintf('Computing superpixels contribution..\n');
    tic
    sp_max_val = max(max(superpixels));
    [~, disp] = min(cv,[],3);
    % apply mask
    disp_mat = disp.*double(matte);
    gimg = rgb2gray(uint8(img));
    num_of_imgs = size(cv, 3);
    cv_sp = zeros(size(cv,1), size(cv,2), size(cv,3));
    THRESH_PEAK2 = 0.4; % indicates minimum ratio to accept second peak
    THRESH_PEAK3 = 0.2; % minimum ratio for third peak
    for i = 1: sp_max_val
        % get pixels
        indices = superpixels == i;
        img = disp_mat(superpixels==i);
        if sum(sum(indices.*matte)) > 0
            %subplot(221)
            %imagesc(double(indices).*double(gimg).*double(matte));
            %axis off
            % divide bins into half of the possibilities
            % if we have 0 - 10 disparities/plane, we want 5 bins
            bins = linspace(0, num_of_imgs, (num_of_imgs));
            %subplot(222)
            %histogram(img, bins);
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
                        ind1 = edges(j);
                    elseif counts(j) == peak2
                        ind2 = edges(j);
                    elseif counts(j) == peak3
                        ind3 = edges(j);
                    end
                end
                x = linspace(0, num_of_imgs, num_of_imgs+1);

                z = z - max(0, 1-((ind1-x)./sigma).^2);
                if peak2 > peak1*THRESH_PEAK2 && size(counts,2)>1
                    z = z - max(0, peak2/peak1-((ind2-x)./sigma).^2);
                end
                if peak3 > peak1*THRESH_PEAK3 && size(counts,2)>2
                    z = z - max(0, peak3/peak1-((ind3-x)./sigma).^2);
                end
                zz = max(z,0);
            end
            %smooth = z;
            %subplot(2,2,[3 4])
            %plot(zz);
            %a = 3;
        else
            zz = ones(1,num_of_imgs+1);
        end
        for k = 1:size(cv,3)
            cv_sp(:,:,k) = cv_sp(:,:,k) + indices.*zz(k);
        end
    end

    fprintf('Done!                       ');
    toc
end