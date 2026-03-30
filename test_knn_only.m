function test_knn_only
inpath = getenv('path2tSpaceInput');
X = csvread(inpath,1,0);
X = double(X);

k = 30;

fprintf('Input: %d cells x %d features\n', size(X,1), size(X,2));

which spdists_knngraph -all

knn = spdists_knngraph(X, k);

fprintf('knn built\n');
fprintf('size(knn) = %d x %d\n', size(knn,1), size(knn,2));
fprintf('nnz(knn) = %d\n', nnz(knn));

if exist('spdists_undirected','file') == 2
    knn2 = spdists_undirected(knn);
    fprintf('undirected graph built\n');
    fprintf('nnz(knn2) = %d\n', nnz(knn2));
end
end