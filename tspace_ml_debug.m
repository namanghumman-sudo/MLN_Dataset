function [allData, tspacem, megaMat, tExplain, pExplain] = tspace_ml_debug(k, l, graph, num_landmarks, tspace_tsne, numPop, perplexity)
% Minimal smoke test for tSpace pipeline
% Purpose: verify CSV input, PCA, kmeans, and CSV output all work
% This intentionally skips spdists_knngraph and runpathFinder

setenv('path2tSpaceInput','C:\Users\ghumm\Downloads\cyt3-master\cyt3-master\src\wanderlust\seurat_export.csv')

setenv('path2tSpaceOutput','C:\Users\ghumm\Downloads\cyt3-master\cyt3-master\src\wanderlust\tspace_allData.csv')

setenv('path2tSpaceOutput2','C:\Users\ghumm\Downloads\cyt3-master\cyt3-master\src\wanderlust\tspace_tExplain.csv')

setenv('path2tSpaceOutput3','C:\Users\ghumm\Downloads\cyt3-master\cyt3-master\src\wanderlust\tspace_pExplain.csv')

getenv('path2tSpaceInput')
getenv('path2tSpaceOutput')

% Read input path from env var
inpath = getenv('path2tSpaceInput');
if isempty(inpath)
    error('path2tSpaceInput not set');
end

out1 = getenv('path2tSpaceOutput');
out2 = getenv('path2tSpaceOutput2');
out3 = getenv('path2tSpaceOutput3');

if isempty(out1) || isempty(out2) || isempty(out3)
    error('Output env vars not set');
end

% Load numeric matrix
X = read_csv_numeric(inpath);
X = double(X);

fprintf('Loaded input: %d cells x %d features\n', size(X,1), size(X,2));

% Basic checks
if size(X,1) < 2 || size(X,2) < 2
    error('Input matrix too small');
end

if nargin < 6 || isempty(numPop)
    numPop = min(5, size(X,1));
end

numPop = min(numPop, size(X,1));

% Kmeans
rng(1);
clusters = kmeans(X, numPop, 'MaxIter', 100, 'Replicates', 1);
fprintf('kmeans finished\n');

% PCA on original data
nPC_p = min(10, size(X,2));
[pScore, pExplained] = pca_svd(X, nPC_p);
pExplain = round(pExplained(:), 2);

% Fake "tspace" matrix for smoke test:
% just reuse first few PCA dimensions
nPC_t = min(5, size(pScore,2));
tspacem = pScore(:,1:nPC_t);

% PCA on fake tspace
[tScore, tExplained] = pca_svd(tspacem, min(3, size(tspacem,2)));
tExplain = round(tExplained(:), 2);

% Build outputs in roughly expected format
index = (1:size(X,1))';
allData = [index, pScore, clusters, tScore, tspacem];
megaMat = [index, clusters, tScore];

% Write outputs
csvwrite(out1, allData);
csvwrite(out2, tExplain);
csvwrite(out3, pExplain);

fprintf('Wrote:\n%s\n%s\n%s\n', out1, out2, out3);
fprintf('DEBUG version completed successfully.\n');

end


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


function [score, explained] = pca_svd(X, nComp)
    X = double(X);
    mu = mean(X, 1);
    Xc = bsxfun(@minus, X, mu);

    [U,S,~] = svd(Xc, 'econ');
    nComp = min(nComp, size(U,2));

    score = U(:,1:nComp) * S(1:nComp,1:nComp);

    latent = diag(S).^2 / max(1, (size(Xc,1)-1));
    explained_all = 100 * latent / sum(latent);
    explained = explained_all(1:nComp);
end