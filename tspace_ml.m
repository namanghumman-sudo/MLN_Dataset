function [allData, tspacem, megaMat, tExplain, pExplain] = tspace_ml(k, l, graph, num_landmarks, tspace_tsne, numPop, perplexity)

if nargin < 1 || isempty(k),             k = 16; end
if nargin < 2 || isempty(l),             l = k; end
if nargin < 3 || isempty(graph),         graph = 1; end
if nargin < 4 || isempty(num_landmarks), num_landmarks = 20; end
if nargin < 5 || isempty(tspace_tsne),   tspace_tsne = 0; end
if nargin < 6 || isempty(numPop),        numPop = 10; end
if nargin < 7 || isempty(perplexity),    perplexity = 30; end
%   - removes dependency on MATLAB pca() (uses SVD PCA)
%   - uses old-style spdists_knngraph calls (positional args only)
%   - robust CSV read whether header exists or not
%
% Env vars expected:
%   path2tSpaceInput   : CSV input path
%   path2tSpaceOutput  : output allData CSV
%   path2tSpaceOutput2 : output tExplain CSV
%   path2tSpaceOutput3 : output pExplain CSV
% Optional:
%   tsp_distMetric     : e.g. euclidean (default)
%   tsp_voting_scheme  : not used here but kept

% --------------------------
% Parameters structure
% --------------------------
gates = cell(2);
parameters.partial_order = [];
parameters.voting_scheme = getenv('tsp_voting_scheme'); %#ok<NASGU>
parameters.band_sample = 1;
parameters.flock_landmarks = 2;
parameters.verbose = 1;            % keep 1 while debugging
parameters.deblur = 0;
parameters.snn = 1;
parameters.search_connected_components = 1;
parameters.knn = [];
parameters.exclude_points = [];
parameters.gates = gates;
parameters.gate_index = 1;
parameters.gateContext = gates{1, 2};
parameters.perplex = perplexity;

parameters.metric = getenv('tsp_distMetric');
if isempty(parameters.metric)
    parameters.metric = 'euclidean';
end

parameters.k = k;
parameters.l = l;
parameters.num_landmarks = num_landmarks;

% --------------------------
% Load input CSV robustly
% --------------------------
inpath = getenv('path2tSpaceInput');
if isempty(inpath)
    error('Environment variable path2tSpaceInput is not set.');
end

sessionData = read_csv_numeric(inpath);

fprintf('tspace_ml: nCells=%d nFeat=%d\n', size(sessionData,1), size(sessionData,2));

% --------------------------
% KMeans populations
% --------------------------
rng(1);
clusters_trajectories = kmeans(sessionData, numPop, 'MaxIter', 10000);

% --------------------------
% PCA on sessionData (SVD, toolbox-free)
% --------------------------
if size(sessionData,2) > 40
    nPC_p = 20;
else
    nPC_p = max(2, round(size(sessionData,2)/2));
end

[pScore, pExplained] = pca_svd(sessionData, nPC_p);
pExplain = round(pExplained(:), 2);

% --------------------------
% Separate cells by kMeans population
% --------------------------
indexPops = zeros(numPop, numPop);
pop = zeros(numPop, 1);
for ii = 1:size(clusters_trajectories,1)
    pop(clusters_trajectories(ii)) = pop(clusters_trajectories(ii)) + 1;
    indexPops(clusters_trajectories(ii), pop(clusters_trajectories(ii))) = ii;
end

% --------------------------
% Build kNN graph (OLD spdists_knngraph calls)
% --------------------------
fprintf('Building kNN graph: k=%d metric=%s\n', parameters.k, parameters.metric);
tic;
knn = knn_graph_compat(sessionData, parameters.k, parameters.metric);
fprintf('kNN computed: %gs\n', toc);

% --------------------------
% Optional deblur
% --------------------------
if (parameters.deblur)
    [i1, j1, ~] = find(knn);
    for ith = 1:size(sessionData,1)
        neigh = j1(i1 == ith);
        if ~isempty(neigh)
            sessionData(ith,:) = median(sessionData(neigh,:), 1);
        end
    end
    knn = knn_graph_compat(sessionData, parameters.k, parameters.metric);
end

% --------------------------
% Shared Nearest Neighbor pruning
% --------------------------
if (parameters.snn ~= 0)
    [j, i, s] = find(knn);
    nData = size(sessionData,1);
    rem = cell(1, nData);

    usePar = license('test','Distrib_Computing_Toolbox');

    if usePar
        parfor ci = 1:nData
            from = (ci-1)*parameters.k+1;
            to   = ci*parameters.k;
            i_inds = from:to;
            i_neighs = j(i_inds);

            localRem = [];
            for i_ind = i_inds
                i_neigh = j(i_ind);

                from2 = (i_neigh-1)*parameters.k+1;
                to2   = i_neigh*parameters.k;
                j_neighs = j(from2:to2);

                if sum(ismember(i_neighs, j_neighs)) < parameters.snn
                    localRem = [localRem, i_ind]; %#ok<AGROW>
                end
            end
            rem{ci} = localRem;
        end
    else
        for ci = 1:nData
            from = (ci-1)*parameters.k+1;
            to   = ci*parameters.k;
            i_inds = from:to;
            i_neighs = j(i_inds);

            for i_ind = i_inds
                i_neigh = j(i_ind);

                from2 = (i_neigh-1)*parameters.k+1;
                to2   = i_neigh*parameters.k;
                j_neighs = j(from2:to2);

                if sum(ismember(i_neighs, j_neighs)) < parameters.snn
                    rem{ci} = [rem{ci} i_ind]; %#ok<AGROW>
                end
            end
        end
    end

    rem = cell2mat(rem);
    i(rem) = []; j(rem) = []; s(rem) = [];
    knn = sparse(j, i, s);
end

% --------------------------
% Make undirected
% --------------------------
if exist('spdists_undirected', 'file') == 2
    knn = spdists_undirected(knn);
else
    knn = max(knn, knn');
end

parameters.knn = knn;
parameters.search_connected_components = true;

% --------------------------
% Run trajectories; average graphs
% --------------------------
tspacem = zeros(size(sessionData,1), numPop);
graph_panel = cell(graph, 1);

for graph_iter = 1:graph

    if (parameters.k ~= parameters.l)
        if exist('spdists_lknn','file') == 2
            lknn = spdists_lknn(knn, parameters.l, parameters.verbose);
        else
            lknn = knn;
        end
    else
        lknn = knn;
    end

    usePar = license('test','Distrib_Computing_Toolbox');
    if usePar
        parfor ii = 1:numPop
            s0 = indexPops(ii,1);
            tspacem(:,ii) = runpathFinder(sessionData, lknn, parameters, s0);
        end
    else
        for ii = 1:numPop
            s0 = indexPops(ii,1);
            tspacem(:,ii) = runpathFinder(sessionData, lknn, parameters, s0);
        end
    end

    graph_panel{graph_iter} = tspacem;
end

tspacem = cat(3, graph_panel{:});
tspacem = mean(tspacem, 3);

% --------------------------
% PCA on tspacem (SVD)
% --------------------------
if size(tspacem,2) > 40
    nPC_t = 20;
else
    nPC_t = max(2, round(size(tspacem,2)/2));
end

[tScore, tExplained] = pca_svd(tspacem, nPC_t);
tExplain = round(tExplained(:), 2);

% --------------------------
% Build outputs
% --------------------------
index = (1:size(sessionData,1))';
allData = cat(2, index, pScore, clusters_trajectories, tScore, tspacem);
megaMat = cat(2, index, clusters_trajectories, tScore);

% --------------------------
% Write outputs
% --------------------------
out1 = getenv('path2tSpaceOutput');
out2 = getenv('path2tSpaceOutput2');
out3 = getenv('path2tSpaceOutput3');

if isempty(out1) || isempty(out2) || isempty(out3)
    error('Output env vars not set: path2tSpaceOutput, path2tSpaceOutput2, path2tSpaceOutput3');
end

csvwrite(out1, allData);
csvwrite(out2, tExplain);
csvwrite(out3, pExplain);

end


% =========================================================
% Helper: robust CSV read (header or no header)
% =========================================================
function X = read_csv_numeric(inpath)
    % Try no-header first
    try
        X = csvread(inpath, 0, 0);
        if any(isnan(X(:))) || isempty(X)
            error('NaN/empty after csvread(0,0)');
        end
        return;
    catch
        % fall through
    end

    % Try skipping 1 header row
    try
        X = csvread(inpath, 1, 0);
        if any(isnan(X(:))) || isempty(X)
            error('NaN/empty after csvread(1,0)');
        end
        return;
    catch ME
        error('Failed to read numeric CSV from %s. Last error: %s', inpath, ME.message);
    end
end


% =========================================================
% Helper: toolbox-free PCA via SVD
% =========================================================
function [score, explained] = pca_svd(X, nComp)
    X = double(X);
    mu = mean(X, 1);
    Xc = bsxfun(@minus, X, mu);

    [U,S,~] = svd(Xc, 'econ');
    nComp = min([nComp, size(U,2)]);

    score = U(:,1:nComp) * S(1:nComp,1:nComp);

    latent = diag(S).^2 / (size(Xc,1)-1);
    explained_all = 100 * latent / sum(latent);
    explained = explained_all(1:nComp);
end


% =========================================================
% Helper: kNN graph builder compatible with older signatures
% =========================================================
function knn = knn_graph_compat(X, k, metric)
    % Try the most common old signatures in order.
    try
        knn = spdists_knngraph(X, k);
        return;
    catch
    end

    try
        knn = spdists_knngraph(X, k, metric);
        return;
    catch ME
        error('spdists_knngraph failed for (X,k) and (X,k,metric). Last error: %s', ME.message);
    end
end