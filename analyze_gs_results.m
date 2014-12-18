% Analyze results of run_sweep_params
rootdir = 'gs_08-Dec-2014 22:26:02';
params_path = fullfile(rootdir, 'params.mat');
results_path = fullfile(rootdir, 'gs_results.mat');

load(params_path);
load(results_path);

%% Determine best parameter setting
accs = [];
for comb_i=1:length(gs_results)
    accs_total = gs_results{comb_i}{4};
    accs = [accs accs_total(1)];
end
[acc_best, ind_best] = max(accs);
accs_total_best = gs_results{ind_best}{4};
fprintf('Best comb_i is: %d (acc: %.4f nb_pos: %d nb_neg: %d)\n', ind_best, acc_best, accs_total_best(2), accs_total_best(3));
for i=1:length(gs_params)
    fprintf('%s: %f\n', gs_params{i}{1}, gs_combs(ind_best,i));
end

%% Visualize
%[labels_pred, params_pred, intermedss, accs_total, accs_digs] = deal(gs_results{4});
labels_pred = gs_results{4}{1};
params_pred = gs_results{4}{2};
intermedss = gs_results{4}{3};
accs_total = gs_results{4}{4};
accs_digs = gs_results{5}{5};
