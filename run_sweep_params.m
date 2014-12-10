% Script to sweep over parameter values
datestr_now = datestr(now);
rootdir = sprintf('gs_%s', datestr_now);
mkdir(rootdir);
partials_dir = fullfile(rootdir, 'partials');
mkdir(partials_dir);
params_path = fullfile(rootdir, 'params.mat');
diary_path = fullfile(rootdir, 'diary.txt');
diary OFF; % make sure any previous diary session is off
diary(diary_path);
fprintf('(do_gridsearch.m) Starting new gridsearch at: %s\n', rootdir);
datestr_now

%% Define parameter combinations
addpath('allcomb');
FLAG_GS = 1; % Signal we're doing a grid search (for do_classify_digit.m)

ANNEAL_SCHED = [...
        %[int NB_BEADS, int MAX_INNER_ITERS, float VAR]
        [7, 4, 0.4]; ...
        [14, 5, 0.2]; ...
        [14, 5, 0.1]; ...
        [14, 16, 0.08]; ...
        [28, 24, 0.04]; ...
        ];
ANNEAL_SCHED = [...
    [8, 5, 0.04]; ...
    [15, 3, 0.02]; ...
    [15, 3, 0.01]; ...
    [15, 10, 0.008]; ...
    [20, 14, 0.0025]; ...
    [45, 25, 0.0006]; ...
    %[15, 28, 0.008]; ...
    %[20, 28, 0.0025]; ...
    %[45, 35, 0.0006]; ...
    [60, 2, 0.0006]; ...
    ];
    
VERBOSE = 0;    
    
gs_params = {...
    {'N_0', [40, 80, 160, 360]}; % N_0 = 160; % "Standard" nb. of pixels. [See pg. 31 of thesis]
    {'PI_N', [0.1, 0.3]}; % PI_N = 0.3; % mixing coef. btwn uniform-noise and gaussian-mixture
    {'C', [1.0]}; % C = 10.0; % weight of E_def. Note: DON'T USE
    };

gs_names = cellfun(@(c) c{1}, gs_params, 'UniformOutput', 0);
gs_vals = cellfun(@(c) c{2}, gs_params, 'UniformOutput', 0);
gs_combs = allcomb(gs_vals{:});

save(params_path, 'gs_params', 'gs_names', 'gs_vals', 'gs_combs', 'datestr_now', 'ANNEAL_SCHED');

%% Do gridsearch
gs_results = {};
tic_gs_total = tic;
for comb_i=1:size(gs_combs, 1)
    fprintf('==== GridSearch (%d/%d) ====\n', comb_i, size(gs_combs,1));
    tic_gs_inner = tic;
    %% Set param values
    curvals = gs_combs(comb_i, :);
    for j=1:length(curvals)
        curname = gs_names{j};
        val = curvals(j);
        eval(sprintf('%s = %d;', curname, val));
        fprintf('%s is: %d\n', curname, val);
    end
    %% Run experiemnt with assigned parameter values
    do_classify_digit;
    %% Save results
    gs_results{comb_i} = {...
        labels_pred, params_pred, intermedss, accs_total, accs_digs, ...
        };
    toc_gs_inner = toc(tic_gs_inner);
    fprintf('Finished comb_i=%d (%.4fs)\n', comb_i, toc_gs_inner);
end
toc_gs_total = toc(tic_gs_total);
fprintf('Finished gridsearch (%.4fs)\n', toc_gs_total);

outpath = fullfile(rootdir, 'gs_results.mat');
save(outpath, 'gs_results');
fprintf('Saved results to: %s\n', outpath);
datestr(now);
diary off;
