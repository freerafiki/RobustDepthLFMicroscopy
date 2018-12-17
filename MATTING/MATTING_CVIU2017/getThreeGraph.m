function L1 = getThreeGraph(I,trimap,K0,K1,K2,K3)
% get three layer graph

[m,n,z] = size(I); N = m*n;

% Computing W1
x = repmat(1:n,m,1);    y = repmat((1:m)',1,n);
featherArray = cat(3,I,y/12,x/12);
X = reshape(featherArray,N,[]);
K = K0;
kdtree = vl_kdtreebuild(X');
[IDX,D] = vl_kdtreequery(kdtree,X',X','NumNeighbors',K+1);
IDX = IDX';		D = D';
IDX = IDX(:,2:end); D = D(:,2:end);

sumD = getsumD(D,IDX,K);
X = [X(:,1:3) X(:,4:5) sumD];
K = K1;

kdtree = vl_kdtreebuild(X');
[IDX,D] = vl_kdtreequery(kdtree,X',X','NumNeighbors',K+1);
IDX = IDX';		D = D';
IDX = IDX(:,2:end); D = D(:,2:end);
[~,IX] = sort(sum(D.^2,2));   IX = IX(1:floor(0.9*N));
flag = zeros(m,n);    flag(IX) = 1;							
neighborsArray = reshape(IDX,m,n,K);
W1 = getGraph(I,trimap,neighborsArray,flag);

% Computing W2 
sumD = getsumD(D,IDX,K);
X = [X(:,1:3) X(:,4:5)./120 sumD];
K = K2;

kdtree = vl_kdtreebuild(X');
[IDX,D] = vl_kdtreequery(kdtree,X',X','NumNeighbors',K+1);
IDX = IDX';		D = D';
IDX = IDX(:,2:end); D = D(:,2:end);
flag = selectPix(I,D,0.15); % originally 0.15
neighborsArray = reshape(IDX,m,n,K);
W2 = getGraph(I,trimap,neighborsArray,flag);

% Computing W3
sumD = getsumD(D,IDX,K);  
X = [X(:,1:3) X(:,4:5)./12 sumD];
K = K3;

kdtree = vl_kdtreebuild(X');
[IDX,D] = vl_kdtreequery(kdtree,X',X','NumNeighbors',K+1);
IDX = IDX';		D = D';
IDX = IDX(:,2:end); D = D(:,2:end);
flag = selectPix(I,D,0.1);
neighborsArray = reshape(IDX,m,n,K);
W3 = getGraph(I,trimap,neighborsArray,flag);

L1 = 0.1*(W3'*W3) + 0.3*(W2'*W2) + W1'*W1;

end
