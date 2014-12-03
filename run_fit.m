trainingImagePath = 'mnist/train-images.idx3-ubyte';
trainingLabelPath = 'mnist/train-labels.idx1-ubyte';
numToTrain = 20;
offset = 0;

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);
I = squeeze(imgs(:,:,6));

level = graythresh(I);

%% Preprocess image
Ip = I;
Ip(Ip < level) = 0;
Ip(Ip ~= 0) = 1;
Ip = bwmorph(Ip,'skel',Inf);

% a silly initial spline just to run through fit_model (simple '2')
% origin to be: [1, 5]
ps = [[1, 5];
    [4, 4.5];
    [3.5, 4.0];
    [2.5, 3.0];
    [1.0, 1.0];
    [2.0, 1.5];
    [4.0, 1.75]];
ps(:, 2) = -ps(:, 2) + ones(size(ps,1), 1) * 6;

% duplicate start/end ctrl pts
ps = [ps(1,:);
    ps;
    ps(end,:)];

tic1 = tic;
[xs, A, t, E_def, E_fit, intermeds] = fit_model(Ip, ps');
toc1 = toc(tic1);
fprintf('Finished fit_model (%.2fs)\n', toc1);

%% Visualize final output
visualize_model(Ip, xs, ps', A, t, 10);
suptitle(sprintf('E_tot: %.2f E_def: %.2f E_fit: %.2f', ...
    E_def+E_fit, E_def, E_fit));
%% Animate iterative process
figure;
for i=1:length(intermeds)
    xs_i = intermeds{i}{1};
    A_i = intermeds{i}{2};
    t_i = intermeds{i}{3};
    E_def_i = intermeds{i}{4};
    E_fit_i = intermeds{i}{5};
    N_B_i = intermeds{i}{6};
    visualize_model(Ip, xs_i, ps', A_i, t_i, N_B_i);
    suptitle(sprintf('Iter %d/%d (E_tot: %.2f E_def: %.2f E_fit: %.2f)', ...
        i, length(intermeds), ...
        E_def_i+E_fit_i, E_def_i, E_fit_i));
    pause
end
