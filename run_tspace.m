
clear functions;
rehash toolboxcache;
rehash path;

% run_tspace.m
here = fileparts(mfilename('fullpath'));
cd(here);

inputCsv = fullfile(here, "seurat_export.csv");
outAll   = fullfile(here, "tspace_allData.csv");
outTE    = fullfile(here, "tspace_tExplain.csv");
outPE    = fullfile(here, "tspace_pExplain.csv");

% parameters
k = 30;          % kNN
l = 15;          % lNN (<=k)
numPop = 10;     % number of start nodes / trajectories
graph_reps = 1;  % keep 1 for now

disp("=== run_tspace started ===");
disp(["Input: ", string(inputCsv)]);
disp(["which(tspace_ml) = ", string(which('tspace_ml'))]);

if ~isfile(inputCsv)
    error("Input CSV not found: %s", inputCsv);
end

[allData, tspacem, megaMat, tExplain, pExplain] = tspace_ml( ...
    inputCsv, outAll, outTE, outPE, k, l, numPop, graph_reps);

disp("Wrote:");
disp(outAll);
disp(outTE);
disp(outPE);
disp("=== run_tspace done ===");