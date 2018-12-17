function feature = getsumD(D,IDX,K)
% get sumD feature

% paramaters

% get sumD feature
layerNum = 120;
sumD = sum(D,2);
for k=1:1:layerNum
	sumD = reshape(D',[],1) + sumD(reshape(IDX',[],1),:);
	sumD = reshape(sumD',K,[]);
	sumD = sumD';
	sumD = sum(sumD,2);
end
sumD = (sumD-min(sumD))./(max(sumD)-min(sumD));

feature = sumD*255/std(sumD);

end
