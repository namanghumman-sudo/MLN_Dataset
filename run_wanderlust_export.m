% run_wanderlust_export.m
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));

disp('Starting run_wanderlust_export');

in_csv   = fullfile(pwd, 'seurat_export.csv');
cell_txt = fullfile(pwd, 'cell_ids.txt');
out_csv  = fullfile(pwd, 'wanderlust_pseudotime.csv');

disp('--- FILE PATHS ---');
disp(in_csv);
disp(cell_txt);
disp(out_csv);

if ~exist(in_csv, 'file')
    error('Input CSV not found: %s', in_csv);
end

if ~exist(cell_txt, 'file')
    error('Cell ID file not found: %s', cell_txt);
end

data = readmatrix(in_csv);

cell_ids = readlines(cell_txt);
cell_ids = cell_ids(cell_ids ~= "");
cell_ids = strtrim(cell_ids);
cell_ids = cell_ids(:);

disp('--- DATA SIZE ---');
disp(size(data));

disp('--- CELL IDS SIZE ---');
disp(size(cell_ids));
disp(numel(cell_ids));

if size(data, 1) ~= numel(cell_ids)
    error('Mismatch: %d rows in data but %d cell IDs.', size(data,1), numel(cell_ids));
end

Options = struct();
Options.metric = 'euclidean';
Options.k = 10;
Options.l = 10;
Options.num_graphs = 1;
Options.s = 1;
Options.num_landmarks = 100;
Options.verbose = 1;
Options.branch = 0;
Options.partial_order = [];
Options.deblur = 0;
Options.snn = 0;
Options.ann = 0;
Options.voting_scheme = 'exponential';
Options.band_sample = 1;
Options.flock_landmarks = 2;
Options.search_connected_components = 1;
Options.plot_landmark_paths = 0;
Options.plot_data = data(:, 1:min(2, size(data,2)));
Options.lnn = [];
Options.landmarks = [];
Options.disallow = [];
Options.cell_clusters = [];
Options.end_clusters = [];
Options.plot_debug_branch = 0;
Options.kEigs = 6;

disp('--- OPTIONS ---');
disp(Options);

G = wanderlust(data, Options);

disp('Reached post-wanderlust section');
disp('--- G fields ---');
g_fields = fieldnames(G);
disp(g_fields);

fprintf('\n--- Top-level field summary ---\n');
for i = 1:numel(g_fields)
    f = g_fields{i};
    val = G.(f);
    if isnumeric(val)
        fprintf('%s -> numeric size %s\n', f, mat2str(size(val)));
    elseif islogical(val)
        fprintf('%s -> logical size %s\n', f, mat2str(size(val)));
    elseif ischar(val)
        fprintf('%s -> char size %s\n', f, mat2str(size(val)));
    elseif isstring(val)
        fprintf('%s -> string size %s\n', f, mat2str(size(val)));
    elseif iscell(val)
        fprintf('%s -> cell size %s\n', f, mat2str(size(val)));
    elseif isstruct(val)
        fprintf('%s -> struct size %s\n', f, mat2str(size(val)));
        subf = fieldnames(val);
        for j = 1:numel(subf)
            sf = subf{j};
            sval = val.(sf);
            if isnumeric(sval)
                fprintf('    %s.%s -> numeric size %s\n', f, sf, mat2str(size(sval)));
            elseif islogical(sval)
                fprintf('    %s.%s -> logical size %s\n', f, sf, mat2str(size(sval)));
            elseif ischar(sval)
                fprintf('    %s.%s -> char size %s\n', f, sf, mat2str(size(sval)));
            elseif isstring(sval)
                fprintf('    %s.%s -> string size %s\n', f, sf, mat2str(size(sval)));
            elseif iscell(sval)
                fprintf('    %s.%s -> cell size %s\n', f, sf, mat2str(size(sval)));
            elseif isstruct(sval)
                fprintf('    %s.%s -> struct size %s\n', f, sf, mat2str(size(sval)));
            else
                fprintf('    %s.%s -> class %s\n', f, sf, class(sval));
            end
        end
    else
        fprintf('%s -> class %s\n', f, class(val));
    end
end

candidate_names = { ...
    'traj', 'trajectory', 'order', 'ordering', 'pseudotime', ...
    'dist', 'dists', 'score', 'scores', 'votes', 'path', 'paths'};

pseudotime = [];
picked_name = '';

for i = 1:numel(candidate_names)
    nm = candidate_names{i};
    if isfield(G, nm)
        val = G.(nm);
        if isnumeric(val) && numel(val) == size(data,1)
            pseudotime = val(:);
            picked_name = ['G.' nm];
            break;
        end
    end
end

if isempty(pseudotime)
    for i = 1:numel(g_fields)
        f = g_fields{i};
        val = G.(f);
        if isnumeric(val) && numel(val) == size(data,1)
            pseudotime = val(:);
            picked_name = ['G.' f];
            break;
        end
    end
end

if isempty(pseudotime)
    for i = 1:numel(g_fields)
        f = g_fields{i};
        val = G.(f);
        if isstruct(val)
            subf = fieldnames(val);
            for j = 1:numel(subf)
                sf = subf{j};
                sval = val.(sf);
                if isnumeric(sval) && numel(sval) == size(data,1)
                    pseudotime = sval(:);
                    picked_name = ['G.' f '.' sf];
                    break;
                end
            end
        end
        if ~isempty(pseudotime)
            break;
        end
    end
end

if isempty(pseudotime)
    error(['Could not automatically identify a pseudotime vector of length ' ...
           num2str(size(data,1)) ...
           '. Inspect the field summary above and set pseudotime manually.']);
end

fprintf('\nSelected pseudotime source: %s\n', picked_name);
fprintf('Pseudotime size: %s\n', mat2str(size(pseudotime)));

if numel(pseudotime) ~= numel(cell_ids)
    error('Mismatch: pseudotime length %d but cell IDs length %d.', numel(pseudotime), numel(cell_ids));
end

disp('Writing output to:');
disp(out_csv);

T = table(cell_ids(:), pseudotime(:), 'VariableNames', {'cell', 'pseudotime'});
writetable(T, out_csv);

disp('File exists after write?');
disp(exist(out_csv, 'file'));

if exist(out_csv, 'file') ~= 2
    error('Output CSV was not created: %s', out_csv);
end

disp('Wanderlust export complete.');
disp(out_csv);

target_dir = 'C:/Users/ghumm/Downloads/cyt3-master/cyt3-master/src/wanderlust';

% force absolute path
out_csv = fullfile(target_dir, 'wanderlust_pseudotime.csv');

disp('--- FINAL WRITE TARGET ---');
disp(out_csv);

% write file
T = table(cell_ids(:), pseudotime(:), 'VariableNames', {'cell', 'pseudotime'});
writetable(T, out_csv);

% verify EXACT location
disp('--- VERIFY WITH FULL PATH ---');
disp(exist(out_csv, 'file'));

% list EXACT directory we wrote to
disp('--- LIST TARGET DIR ---');
disp(dir(target_dir));

fid = fopen(fullfile(target_dir, 'matlab_test_write.txt'), 'w');
fprintf(fid, 'hello\n');
fclose(fid);

disp(exist(fullfile(target_dir, 'matlab_test_write.txt'), 'file'))
system(['dir /a "', target_dir, '"'])