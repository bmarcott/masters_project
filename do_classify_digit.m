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
    LAMBDA_REG = 0.0; % regularization param for min_E_def()
    P_A = 5.0; % for classification: E_tot = E_def + P_A*E_fit
    VERBOSE = 1;
    RNG_SEED = 42;
end
if 1
    tmp_fpath = sprintf('%s_diary.txt', datestr(now, 'dd-mmm-yyyy_HH_MM_SS'));
    diary(tmp_fpath);
end
%% Read Images
rng(RNG_SEED);
trainingImagePath = 'mnist/train-images.idx3-ubyte';
trainingLabelPath = 'mnist/train-labels.idx1-ubyte';
numToTrain = 10000;
offset = 0;

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);
%% Decide which images to use
inds_2 = find(labels==2)';
inds_3 = find(labels==3)';
NB_EXS = 200;
inds_2 = inds_2(randperm(length(inds_2), NB_EXS/2));
inds_3 = inds_3(randperm(length(inds_3), NB_EXS/2));
inds_imgs = [inds_2, inds_3];

%inds_imgs = [199];
%inds_imgs = inds_imgs(10);

% A few 'easy' examples
%inds_2 = [6, 17, 26, 29, 77, 83, 110, 118, 121, 123, 144];
%inds_3 = [8, 11, 13, 28, 45, 50, 51, 75, 87];
%inds_imgs = [inds_2, inds_3];

%inds_imgs = [26]; % an easy 2
%inds_imgs = [77]; % an easy 2
%% Visualize the images
if 0
    hfig = figure;
    for ind=inds_imgs
        I = imresize(squeeze(imgs(:,:,ind)), 5.0);
        imshow(I);
        title(sprintf('Img ID: %d', ind));
        pause;
    end
    close(hfig);
end

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

% curly two
twoCPs = [[1.5, 3.5];
    [2.4, 3.9];
    [2.4, 2.75]
    [1.5, 1.5];
    [1.5, 2.5];
    [1.85, 2];
    [2, 1.9];
    [2.4, 1.5]];
twoCPs(:, 2) = -twoCPs(:, 2) + ones(size(twoCPs,1), 1) * 6;
% duplicate start/end ctrl pts
twoCPs = [twoCPs(1,:);
    twoCPs;
    twoCPs(end,:)];
% Digit model 3
threeCPs = [[1.25, 3.55];
      [1.6, 4.25];
      [1.7, 2.65];
      [1.25, 2];
      [1.5, 2.5];
      [1.65, 1.1];
      [1.6, 0];
      [1.2, 0.5]];
threeCPs(:, 2) = -threeCPs(:, 2) + ones(size(threeCPs,1), 1) * 6;  
threeCPs = [threeCPs(1,:);
            threeCPs;
            threeCPs(end,:)];
models = {...
    {2, twoCPs'}, ...
    {3, threeCPs'}, ...
    };

if (~exist('FLAG_GS', 'var'))
    %% Visualize spline models
    hfig = figure;
    for i=1:length(models)
        % Compute affine-trans that maps spline-ctrl-pts to fill up image
        ps = models{i}{2}; % [2 x N]
        x1 = min(ps(1,:));
        y1 = min(ps(2,:));
        x2 = max(ps(1,:));
        y2 = max(ps(2,:));
        w = x2 - x1; h = y2 - y1;
        bbox_obj = [[x1 y1];
            [x1+w, y1];
            [x1, y1+h];
            [x1+w, y1+h]];
        bbox_img = [[0 0];
            [20, 0];
            [0, 20];
            [20, 20]];
        tform = cp2tform(bbox_obj, bbox_img, 'affine');
        % Note: T is arranged in a weird way (transpose it to get familiar affine)
        T = tform.tdata.T'; % transpose
        tmpA = T(1:2, 1:2);
        tmpt = T(1:2, 3);
        ps_ = tmpA*ps + repmat(tmpt, [1, size(ps, 2)]);
        subplot(1, length(models), i);
        Ishow = rasterize_spline(ps_, [20,20]);
        imshow(Ishow);
        title(sprintf('"%d" model', models{i}{1}));
    end
end

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

%% Save workspace to file.
if 1
tmp_fpath = sprintf('%s.mat', datestr(now, 'dd-mmm-yyyy_HH_MM_SS'));
clear('imgs');
save(tmp_fpath);
fprintf('Saved results to: %s\n', tmp_fpath);
end

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

%% Do alternate classification
labels_pred2 = [];
for i=1:length(intermedss) % for each image
    label_pred = labels_pred(i);
    label_true = labels(inds_imgs(i));
    Es = [];
    for m_i=1:length(intermedss{i})
        % {xs_est, A, t, E_def, E_fit, N_B}
        res = intermedss{i}{m_i}{end};
        Es(m_i) = res{4} + P_A*res{5};
    end
    [~, ind_pred] = min(Es);
    labels_pred2(i) = models{ind_pred}{1};
end
nb_pos2 = sum(labels(inds_imgs) == labels_pred2');
nb_neg2 = length(inds_imgs) - nb_pos2;
acc2 = nb_pos2 / (nb_pos2 + nb_neg2);
accs_total2 = [acc2, nb_pos2, nb_neg2, P_A];
fprintf('Alternate Classification (P_A=%f):\n', P_A);
fprintf('    Acc: %.4f (nb_pos: %d nb_neg: %d)\n', acc2, nb_pos2, nb_neg2);

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
%% Visualize iterations
for i = 1 : length(params_pred) % for each image
    % Visualize final output
    xs = params_pred{i}{1};
    A = params_pred{i}{2};
    t = params_pred{i}{3};
    E_def = params_pred{i}{4};
    E_fit = params_pred{i}{5};
    intermeds = intermedss{i};
    Ip = imgs_p(:,:,i);
%     visualize_model(Ip, xs, ps', A, t, 10);
    suptitle(sprintf('E_tot: %.2f E_def: %.2f E_fit: %.2f', ...
    E_def+E_fit, E_def, E_fit));
    % Animate iterative process
    h = figure;
    for model_ind=1:length(intermeds) % for each digit/model
        ps = models{model_ind}{2}';
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
