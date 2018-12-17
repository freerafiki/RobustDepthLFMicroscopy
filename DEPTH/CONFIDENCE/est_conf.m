function [confidence, ind] = est_conf(cost_volume)

confidence = zeros(size(cost_volume,1), size(cost_volume,2));
[cost,ind] = min(cost_volume,[],3);
d_max = size(cost_volume, 3);
d_min = 1;
cost = squeeze(cost);
ind = squeeze(ind);

%cd DEPTH/CPP/
%mex calc_conf.cpp
%cd ../..
%conf = calc_conf(cost_volume, ind);
tic
f1 = (d_max - d_min) / 3.0;
sum = 0.0;
for i=1:size(confidence, 1)-1
    for j=1:size(confidence, 2)-1
        sum = 0.0;
        for d=d_min:d_max-1
            delta_d = abs(d - ind(i,j));
            delta_c = cost_volume(i,j,d) - cost(i,j);
            c_mean = mean(cost_volume(i,j,:));
            sum = sum + ((max(min(delta_d - 1, f1),0)^2)/(max(delta_c-c_mean/5,1))); 
        end
        confidence(i,j) = 1 / sum;
    end
    %fprintf('Analyzing %d-th row, %d%% completed.. \n',i,round(i/size(confidence, 1)*100));
end
toc

end
