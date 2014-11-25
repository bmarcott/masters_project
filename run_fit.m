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
ps = [[0, 0];
    [1, 2];
    [1.5, 1.75];
    [2.5, 0.1];
    [3.5, 3.0];
    [4.5, 4.0];
    [5.5, 1.0];
    [6.5, 0.0];
    ];

% duplicate start/end ctrl pts
ps = [ps(1,:);
    ps;
    ps(end,:)];

[xs, A, t, E_def, E_fit] = fit_model(Ip, ps');
