function W = getGraph(I,trimap,neighborsArray,flag)
% computing weighing matrix

[m,n,z] = size(I); N = m*n;		K = size(neighborsArray,3);
tol=1e-3; % regularlizer in case constrained fits are ill conditioned

featherArray = cat(3,I);
X = reshape(featherArray,N,[]);
w_array = zeros(K,N);
i_array = ones(K,N);
j_array = ones(K,N);

indArray = find(~(trimap == 255 | trimap == 0 | flag == 0));
neighborsIndex = reshape(neighborsArray,N,K);
neighbors = X(reshape(neighborsIndex',1,[]),:);
neighbors = reshape(neighbors',[size(X,2),K,N]);
A = neighbors(:,1,:) - permute(X,[2,3,1]);					% shift ith pt to origin
w = zeros(K,1);

for ind = indArray'
                     
    C = A(:,:,ind)'*A(:,:,ind);                    % local covariance 
    if trace(C) ==  0
        %continue;
        w = ones(K,1);
    else
        C = C + eye(K,K)*tol*trace(C);               % regularlization (K>D)
        w = C\ones(K,1);                             % solve Cw=1
    end
   
    w = w/sum(w);                                    % enforce sum(w)=1
    w_array(:,ind) = w;
end
i_array(:,indArray) = ones(K,1)*indArray';
j_array(:,indArray) = neighborsIndex(indArray,:)'; 

i_array = i_array(:);
j_array = j_array(:);
w_array = w_array(:);
o_array = flag(:);
O = spdiags(o_array,0,N,N);
W = speye(N) - sparse(i_array,j_array,w_array,N,N);
W = O*W;
end

