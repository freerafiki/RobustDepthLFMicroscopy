function [depth, conf] = extract_depth_and_confidence(cv, c_img, matte, which_depth, which_conf)

if which_depth == 0
    
    [~, depth] = min(cv, [], 3);
    
end

depth = depth .* double(matte);

conf = single_confidence(cv, c_img, depth, which_conf, matte);