function test_runpathFinder
% Simple smoke test for runpathFinder on one seed
% Requires env var path2tSpaceInput set and runpathFinder on path

% ---- config ----
k = 30;          % knn neighbors
useVerbose = true;

setenv('path2tSpaceInput','C:\Users\ghumm\Downloads\cyt3-master\cyt3-master\src\wanderlust\seurat_export.csv')

% ---- checks ----
inpath = getenv('path2tSpaceInput');
if isempty(inpath)
    error('path2tSpaceInput not set (use setenv)');
end

% show where runpathFinder and related functions are
disp('Which runpathFinder:'); which runpathFinder -all
disp('Which spdists_knngraph:'); which spdists_knngraph -all
disp('Which spdists_undirected:'); which spdists_undirected -all

% ---- read data ----
X = read_csv_numeric(inpath);
X = double(X);
nCells = size(X,1);
fprintf('Loaded %d cells x %d features\n', nCells, size(X,2));

% ---- kmeans quick cluster assignment (small numPop) ----
numPop = min(6, max(2, round(nCells/1000))); % small number for testing
rng(1);
clusters = kmeans(X, numPop, 'MaxIter', 100, 'Replicates', 1);
fprintf('kmeans -> %d clusters\n', numPop);

% pick first non-empty cluster seed
s0 = [];
for c = 1:numPop
    idx = find(clusters == c, 1, 'first');
    if ~isempty(idx)
        s0 = idx;
        fprintf('Using seed s0 = %d from cluster %d\n', s0, c);
        break;
    end
end
if isempty(s0)
    error('No non-empty cluster found');
end

% ---- build knn graph ----
fprintf('Building knn (k=%d)\n', k);
knn = spdists_knngraph(X, k);         % old signature tested earlier
fprintf('knn built: size = %d x %d, nnz = %d\n', size(knn,1), size(knn,2), nnz(knn));

% make undirected
if exist('spdists_undirected','file') == 2
    knn = spdists_undirected(knn);
else
    knn = max(knn, knn');
end

% ---- create minimal parameters struct used by runpathFinder ----
parameters = struct();
parameters.k = k;
parameters.l = k;            % keep same for test
parameters.num_landmarks = 10;
parameters.verbose = 1;
parameters.flock_landmarks = 2;
parameters.deblur = 0;
parameters.snn = 1;
parameters.search_connected_components = 1;
parameters.knn = knn;
parameters.metric = 'euclidean';

% ---- call runpathFinder once ----
fprintf('Calling runpathFinder for seed %d ...\n', s0);
try
    pathVec = runpathFinder(X, knn, parameters, s0);
    fprintf('runpathFinder returned vector length: %d\n', numel(pathVec));
    fprintf('Summary: min=%g max=%g mean=%g nnz=%d\n', min(pathVec), max(pathVec), mean(pathVec), nnz(pathVec));
    % show first 10 values
    disp('first 10 values of pathVec:');
    disp(pathVec(1:min(10,end))');
catch ME
    fprintf('runpathFinder failed with error:\n%s\n', ME.message);
    fprintf('Call stack:\n');
    for kf = 1:numel(ME.stack)
        fprintf(' %s line %d\n', ME.stack(kf).name, ME.stack(kf).line);
    end
end

end


% -------------------------
% Helper: robust CSV read
% -------------------------
function X = read_csv_numeric(inpath)
    try
        X = csvread(inpath, 0, 0);
        if isempty(X) || any(isnan(X(:)))
            error('csvread(0,0) produced empty/NaN');
        end
        return;
    catch
    end
    try
        X = csvread(inpath, 1, 0);
        if isempty(X) || any(isnan(X(:)))
            error('csvread(1,0) produced empty/NaN');
        end
        return;
    catch ME
        error('Failed to read CSV: %s', ME.message);
    end
end