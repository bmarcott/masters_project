function [label, params_out, intermedss] = classify_digit(I, models, anneal_sched, N_0, pi_n, C, verbose)
%CLASSIFY_DIGIT
%INPUT
%  matrix I: [h x w]
%  cell models: 
%    Stores models in the following way:
%      models{i} -> {str label, matrix cs}
%    where cs is [2 x n]
%OUTPUT
%  int label:
%  cell params:
%    params{i} -> {matrix xs_est, matrix A, t, E_def, E_fit}
%  cell intermedss:
%    intermedss{model} -> {xs_est, A, t, E_def, E_fit}
Es = [];
params = {};
intermedss = {}; % 
for i=1:length(models)
    if verbose
        lbl = models{i}{1};
        fprintf('[i=%d/%d] Fitting model "%d" to image.\n', i, length(models), lbl);
    end
    cs = models{i}{2};
    [xs_est, A_i, t_i, E_def_i, E_fit_i, intermeds] = fit_model(I, cs, anneal_sched, N_0, pi_n, C, verbose);
    Es(i) = C*E_def_i + E_fit_i;
    params{i} = {xs_est, A_i, t_i, E_def_i, E_fit_i};
    intermedss{i} = intermeds;
end

[min_val, ind_min] = min(Es);
label = models{ind_min}{1};
params_out = params{ind_min};
end
