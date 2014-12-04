%% Read Images
trainingImagePath = 'mnist/train-images.idx3-ubyte';
trainingLabelPath = 'mnist/train-labels.idx1-ubyte';
numToTrain = 100;
offset = 0;

inds_2 = [6, 17, 26, 29, 77, 83];
inds_3 = [8, 11, 13, 28, 45, 50, 51, 75, 87];

inds_2 = inds_2(1:2);
inds_3 = inds_3(1:2);

inds_2 = [6];
inds_3 = [];

inds_imgs = [inds_2, inds_3];

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);

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
threeCPs = [threeCPs(1,:);
            threeCPs;
            threeCPs(end,:)];
models = {{2, twoCPs'}, ...
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
    [labels_pred(i), params_pred{i}, intermedss{i}] = classify_digit(Ip, models, 1);
end
toc1 = toc(tic1);
fprintf('Done fitting. (%.2fs)\n', toc1);

%% Evaluate
nb_pos = sum(labels(inds_imgs) == labels_pred');
nb_neg = length(inds_imgs) - nb_pos;
acc = nb_pos / (nb_pos + nb_neg);
fprintf('Acc: %.4f (nb_pos: %d nb_neg: %d)\n', acc, nb_pos, nb_neg);

fprintf('Done.\n');

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
    %for model_ind=1:length(intermeds) % for each digit/model
    for model_ind=2
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
