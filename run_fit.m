trainingImagePath = 'mnist/train-images.idx3-ubyte';
trainingLabelPath = 'mnist/train-labels.idx1-ubyte';
numToTrain = 20;
offset = 0;

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);
I = squeeze(imgs(:,:,1));

level = graythresh(I);

%% Preprocess image
Ip = I;
Ip(Ip < level) = 0;
Ip(Ip ~= 0) = 1;
Ip = bwmorph(Ip,'skel',Inf);

% a silly initial spline just to run through fit_model
ps = [[1, 1];
    [1.2, 2];
    [1.5, 1.75];
    [2.5, 1.1];
    [3.5, 3.0];
    [4.5, 4.0];
    [5.5, 1.0];
    [6.5, 2.2];
    ];

% duplicate start/end ctrl pts
ps = [ps(1,:);
    ps;
    ps(end,:)];

[xs, A, t, E_def, E_fit, intermeds] = fit_model(Ip, ps');

%% Visualize final output
visualize_model(Ip, xs, ps', A, t);
suptitle(sprintf('E_tot: %.2f E_def: %.2f E_fit: %.2f', ...
    E_def+E_fit, E_def, E_fit));

%% Animate iterative process
for i=1:length(intermeds)
    xs_i = intermeds{i}{1};
    A_i = intermeds{i}{2};
    t_i = intermeds{i}{3};
    E_def_i = intermeds{i}{4};
    E_fit_i = intermeds{i}{5};
    hFig = visualize_model(Ip, xs_i, ps', A_i, t_i);
    suptitle(sprintf('Iter %d/%d (E_tot: %.2f E_def: %.2f E_fit: %.2f)', ...
        i, length(intermeds), ...
        E_def_i+E_fit_i, E_def_i, E_fit_i));
    pause
    close(hFig);
end
