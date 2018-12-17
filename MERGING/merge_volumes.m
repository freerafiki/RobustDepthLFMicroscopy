function merged_cv = merge_volumes(cv1, cv2, phi, conf1, conf2, matte, priorsmap, solution_type)


if solution_type == 0
    
    %% CONFIDENCE ASSESSMENT
    % use confidence to estimate a weight for the two volumes (use the
    % average) not so clever
    sum_confD = sum(sum(conf1.*double(matte))); % sum of values wihtout masked part
    num_of_valsD = sum(sum(matte));
    mean_confD = sum_confD / num_of_valsD;
    sum_confC = sum(sum(conf2.*double(matte))); % sum of values wihtout masked part
    num_of_valsC = sum(sum(matte));
    mean_confC = sum_confC / num_of_valsC;
    %std_confD = std(std(confD));
    %std_confC = std(std(confC));


    %% PICK PHI
    if phi == -1
        % choose phi
        phi = mean_confD / (mean_confD + mean_confC);
    end

    merged_cv = phi.*cv1 + (1-phi).*cv2;

elseif solution_type == 1
    
    %% IDEA 1 
    % here we create the data term for the optimization
    % so E = C_c + C_d + k x |C_c - C_d|
    ro = max(max(max(max(cv1))), max(max(max(cv2))) ) / 10;
    merged_cv = cv1 + cv2 + ro .* abs(cv1 - cv2);

elseif solution_type == 2
    
    %% IDEA 2
    % here we create the data term for the optimization
    % so E = a x C_c + (1-a) x C_d + k x |C_c - C_d|
    % where a is calculated using peak ratio
    sorted = sort(cv1, 3); % MINIMUM VALUES
    pkr1 = min(1,sorted(:,:,1)./ sorted(:,:,2)); % NAIVE PEAK RATIO
    sorted2 = sort(cv2, 3); 
    pkr2 = min(1,sorted2(:,:,1)./sorted2(:,:,2)); % NAIVE PEAK RATIO
    a = pkr2 ./ (pkr1 + pkr2);
    ro = max(max(max(max(cv1))), max(max(max(cv2))) ) / 5;
    merged_cv = zeros(size(cv1));
    for j = 1:size(merged_cv, 3)
        a1 = a .* cv1(:,:,j) + (1 - a) .* cv2(:,:,j);
        ro1 = ro .* abs(cv1(:,:,j) - cv2(:,:,j));  
        slice = a .* cv1(:,:,j) + (1 - a) .* cv2(:,:,j) + ro .* abs(cv1(:,:,j) - cv2(:,:,j));  
        %subplot(131), imagesc(a1), subplot(132), imagesc(ro1), subplot(133), imagesc(slice);
        merged_cv(:,:,j) = a .* cv1(:,:,j) + (1 - a) .* cv2(:,:,j) + ro .* abs(cv1(:,:,j) - cv2(:,:,j));  
    end
    
elseif solution_type == 3
    
    %% IDEA 3
    % here we create the data term for the optimization
    % so E = a x C_c + (1-a) x C_d + k x |C_c - C_d|
    % where a is calculated using the curve around the minimum
    [cost1, min1] = min(cv1, [], 3); % MINIMUM VALUES
    [cost2, min2] = min(cv2, [], 3); % MINIMUM VALUES
    pkr1 = zeros(size(min1, 1), size(min1, 2));
    pkr2 = zeros(size(min1, 1), size(min1, 2));
    for i = 1:size(cv1, 1)
        for j = 1:size(cv1,2)
            pkr1(i,j) = min(cv1(i,j,min1(i,j))./cv1(i,j,min(size(cv1,3), min1(i,j)+1)), cv1(i,j,min1(i,j))./cv1(i,j,max(1, min1(i,j)-1)));
            pkr2(i,j) = min(cv2(i,j,min2(i,j))./cv2(i,j,min(size(cv1,3), min2(i,j)+1)), cv2(i,j,min2(i,j))./cv2(i,j,max(1, min2(i,j)-1)));
        end
    end
    a = pkr1 ./ (pkr1 + pkr2);
    ro = 0; % max(max(max(max(cv1))), max(max(max(cv2))) ) / 5;
    merged_cv = zeros(size(cv1));
    for j = 1:size(merged_cv, 3)
        a1 = a .* cv1(:,:,j) + (1 - a) .* cv2(:,:,j);
        ro1 = ro .* abs(cv1(:,:,j) - cv2(:,:,j));  
        slice = a .* cv1(:,:,j) + (1 - a) .* cv2(:,:,j) + ro .* abs(cv1(:,:,j) - cv2(:,:,j));  
        %subplot(131), imagesc(a1), subplot(132), imagesc(ro1), subplot(133), imagesc(slice);
        merged_cv(:,:,j) = a .* cv1(:,:,j) + (1 - a) .* cv2(:,:,j) + ro .* abs(cv1(:,:,j) - cv2(:,:,j));  
    end
        
elseif solution_type == 4
    
    %% IDEA 4
    % we don't need to calculate the confidence
    % the weight for the two data terms can be generated based on how much
    % they "agree" on a minimum - if the two temptative minimum are very
    % far between them we stay with defocus guess (intuitively this case is
    % the case for less reliable points, where defocus should have at least
    % guessed the right depth slice), while if they are close we stay with
    % correspondence guess (correspondence is usually, if correct, more
    % precise than defocus that will have low values on a larger depth
    % slice). We model this by using weight a that is build by factor x
    % difference between the two indices.
    %
    % Moreover, we want to strengthen the good matches. So, we find the
    % points where both defocus and stereo agrees and we reinforce the
    % minimum there (we increase the costs of the other labels). By doing
    % so, in the optimization step it would be harder to change their
    % label, so they would work as "ground control points".
    
    % temptative minimum
    [defcost, defmin] = min(cv1, [], 3);
    [smcost, smmin] = min(cv2, [], 3);
    a = imabsdiff(defmin, smmin) ./ size(cv1,3);  % a = [0, max_disp]
    gcp = (a == 0); % these are the good points, where they agree
    penalty_factor = 0.1;
    merged_cv = zeros(size(cv1));
    for j = 1:size(merged_cv, 3)
        % at index j, we want to penalize the cost c(j) by a quantity that
        % is abs(correct_index - j), so that if correct_index == j we have
        % the minimum, and the other values are penalized. These should
        % happen only for ground control points (so .* gcp) so we use the
        % mask
        to_be_penalized = imabsdiff((ones(size(gcp)).*j), defmin.*gcp)./ size(cv1,3) .* (gcp);
        merged_cv(:,:,j) = a .* cv1(:,:,j) + (1 - a) .* cv2(:,:,j) + penalty_factor .* to_be_penalized;
    end
end