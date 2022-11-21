function M = metzler_indexer(n)

%find the off-diagonal elements of a matrix

% nt = n^2-n;
M = logical(reshape(ones(n)-eye(n), [], 1));
% M = find(reshape(ones(n)-eye(n), [], 1));
% i = (1:nt);
% v = ones(nt,1);

% M = sparse(i, j, v, nt, n^2);

end