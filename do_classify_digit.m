%% Define parameters
if (~exist('FLAG_GS', 'var'))
    ANNEAL_SCHED = [...
        %[int NB_BEADS, int MAX_INNER_ITERS, float VAR]
        [8, 5, 0.4]; ...
        [16, 12, 0.3]; ...
        [32, 14, 0.2]; ...
        [50, 16, 0.1]; ...
        [50, 20, 0.08]; ...
        [60, 30, 0.08]; ...
        ];
    ANNEAL_SCHED = [...
        [8, 5, 0.04]; ...
        [15, 3, 0.02]; ...
        [15, 3, 0.01]; ...
        [15, 10, 0.008]; ...
        [20, 12, 0.0025]; ...
        [30, 20, 0.0006]; ...
        %[15, 28, 0.008]; ...
        %[20, 28, 0.0025]; ...
        %[45, 35, 0.0006]; ...
        [60, 2, 0.0006]; ...
        ];
    N_0 = 160; % "Standard" nb. of pixels. [See pg. 31 of thesis]
    PI_N = 0.3; % mixing coef. btwn uniform-noise and gaussian-mixture
    C = 1.0; % weight of E_def
    LAMBDA_REG = 10.0; % regularization param for min_E_def()
    VERBOSE = 1;
    RNG_SEED = 42;
end

%% Read Images
rng(RNG_SEED);
trainingImagePath = 'mnist/train-images.idx3-ubyte';
trainingLabelPath = 'mnist/train-labels.idx1-ubyte';
numToTrain = 10000;
offset = 0;

%inds_2 = [6, 17, 26, 29, 77, 83, 110, 118, 121, 123, 144];
%inds_3 = [8, 11, 13, 28, 45, 50, 51, 75, 87];

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);
inds_2 = find(labels==2)';
inds_3 = find(labels==3)';
NB_EXS = 50;
inds_2 = inds_2(randperm(length(inds_2), NB_EXS/2));
inds_3 = inds_3(randperm(length(inds_3), NB_EXS/2));

inds_imgs = [inds_2, inds_3];
%inds_imgs = [199];
inds_imgs = inds_imgs(10);
%% Define Models
% Digit model 2
twoCPs = [[1.5, 3];
    [2.2, 4];
    [3.2, 3.9]
    [4, 3];
    [2.5, 2.5];
    [1.25, 1.5];
    [2.0, 1.5];
    [4.0, 1.5]];
twoCPs(:, 2) = -twoCPs(:, 2) + ones(size(twoCPs,1), 1) * 6;
% duplicate start/end ctrl pts
twoCPs = [twoCPs(1,:);
    twoCPs;
    twoCPs(end,:)];
% Digit model 3
threeCPs = [[1, 4.5];
      [2, 4];
      [1.5, 2.5];
      [1.25, 2.5];
      [1.5, 2.5];
      [2, 1];
      [1, 0.5];
      [1, 0.5]];
threeCPs(:, 2) = -threeCPs(:, 2) + ones(size(threeCPs,1), 1) * 6;  
threeCPs = [threeCPs(1,:);
            threeCPs;
            threeCPs(end,:)];
models = {...
    {2, twoCPs'}, ...
    {3, threeCPs'}, ...
    };

%% Preprocess imgs
imgs_p = [];
for i=1:length(inds_imgs)
    ind = inds_imgs(i);
    I = imgs(:,:,ind);
    imgs_p(:,:,i) = preprocess_img(I);
end
 
%% Do classification
tic1 = tic;
labels_pred = [];
params_pred = {};
intermedss = {};
for i=1:size(imgs_p, 3)
    fprintf('== Classifying image %d/%d ==\n', i, size(imgs_p, 3));
    Ip = squeeze(imgs_p(:,:,i));
    tic_inner = tic;
    [labels_pred(i), params_pred{i}, intermedss{i}] = classify_digit(Ip, models, ANNEAL_SCHED, N_0, PI_N, C, LAMBDA_REG, VERBOSE);
    toc_inner = toc(tic_inner);
    fprintf('== Finished image %d/%d (%.4fs) ==\n', i, size(imgs_p,3), toc_inner);
end
toc1 = toc(tic1);
fprintf('Done fitting. (%.2fs)\n', toc1);

%% Evaluate
nb_pos = sum(labels(inds_imgs) == labels_pred');
nb_neg = length(inds_imgs) - nb_pos;
acc = nb_pos / (nb_pos + nb_neg);
fprintf('Acc: %.4f (nb_pos: %d nb_neg: %d)\n', acc, nb_pos, nb_neg);
accs_total = [acc, nb_pos, nb_neg];

% Digit-level accuracies
accs_digs = [];
for i=1:length(models)
    model = models{i};
    digit_cur = model{1};
    inds_dig = find(labels(inds_imgs) == digit_cur); 
    nb_pos_dig = sum(digit_cur == labels_pred(inds_dig));
    nb_neg_dig = length(inds_dig) - nb_pos_dig;
    acc_dig = nb_pos_dig/(nb_pos_dig+nb_neg_dig);
    accs_digs(i,:) = [acc_dig, nb_pos_dig, nb_neg_dig];
    fprintf('[Digit %d] Acc: %.4f (nb_pos: %d nb_neg: %d)\n', digit_cur, acc_dig, nb_pos_dig, nb_neg_dig);
end
fprintf('Done.\n');

if ~exist('FLAG_GS', 'var') % don't visualize if we're doing gridsearch
%% Visualize classifications
for i=1:length(inds_imgs)
    Ip = squeeze(imgs_p(:,:,i));
    if labels(inds_imgs(i)) == 2
        ps = models{1}{2}';
    else
        ps = models{2}{2}';
    end
    intermeds = intermedss{i}; % intermeds for image i
    hfigs = [];
    for dig_i=1:length(models)
        intermeds_dig = intermeds{dig_i}; % intermeds for digit dig_i
        xs = intermeds_dig{end}{1};
        A = intermeds_dig{end}{2};
        t = intermeds_dig{end}{3};
        E_def = intermeds_dig{end}{4};
        E_fit = intermeds_dig{end}{5};
        digit_cur = models{dig_i}{1};
        hfig = figure;
        hfigs = [hfigs hfig];
        visualize_model(Ip, xs, ps', A, t, ANNEAL_SCHED(end,1)); hold on;
        suptitle(sprintf('Img %d [%d/%d] Model "%d" Pred: %d True: %d E_tot= %.2f*%.2f + %.2f = %.2f', ...
            inds_imgs(i), i, length(inds_imgs), digit_cur, labels_pred(i), labels(inds_imgs(i)), ...
            C, E_def, E_fit, C*E_def+E_fit));
    end
    pause;
    close(hfigs);
end
end

if 0
%% Visualize
for i = 1 : length(params_pred) % for each image
    % Visualize final output
    xs = params_pred{i}{1};
    A = params_pred{i}{2};
    t = params_pred{i}{3};
    E_def = params_pred{i}{4};
    E_fit = params_pred{i}{5};
    intermeds = intermedss{i};
    ps = models{i}{2}';
    Ip = imgs_p(:,:,i);
    visualize_model(Ip, xs, ps', A, t, 10);
    suptitle(sprintf('E_tot: %.2f E_def: %.2f E_fit: %.2f', ...
    E_def+E_fit, E_def, E_fit));
    % Animate iterative process
    figure;
    for model_ind=1:length(intermeds) % for each digit/model
        intermed_mdl = intermeds{model_ind}; 
        for step_i=1:length(intermed_mdl); % for each step/iter
            xs_i = intermed_mdl{step_i}{1};
            A_i = intermed_mdl{step_i}{2};
            t_i = intermed_mdl{step_i}{3};
            E_def_i = intermed_mdl{step_i}{4};
            E_fit_i = intermed_mdl{step_i}{5};
            N_B_i = intermed_mdl{step_i}{6};
            visualize_model(Ip, xs_i, ps', A_i, t_i, N_B_i);
            suptitle(sprintf('Iter %d/%d (E_tot: %.2f E_def: %.2f E_fit: %.2f)', ...
                step_i, length(intermed_mdl), ...
                E_def_i+E_fit_i, E_def_i, E_fit_i));
            pause
        end
    end
end
end
