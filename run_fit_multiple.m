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
ps = [[1.5, 3];
    [2.2, 4];
    [3.2, 3.9]
    [4, 3];
    [2.5, 2.5];
    [1.25, 1.5];
    [2.0, 1.5];
    [4.0, 1.5]];
ps(:, 2) = -ps(:, 2) + ones(size(ps,1), 1) * 6;

% duplicate start/end ctrl pts
ps = [ps(1,:);
    ps;
    ps(end,:)];

% another silly initial spline, a simple '3'
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
        
threeCPs(:, 2) = -threeCPs(:, 2) + ones(size(threeCPs,1), 1) * 6;
labels = ['2'; '3'];
CPs = {ps; threeCPs};
xss = cell(size(labels,1),1);
As = cell(size(labels,1),1);
ts = cell(size(labels,1),1);
E_defs = zeros(size(labels,1),1);
E_fits = zeros(size(labels,1),1);
intermedss = cell(size(labels,1),1);
energies = zeros(size(labels,1),1);

for i = 1 : size(labels,1)
    [xs, A, t, E_def, E_fit, intermeds] = fit_model(Ip, CPs{i}');
    xss{i} = xs;
    As{i} = A;
    ts{i} = t;
    E_defs(i) = E_def;
    E_fits(i) = E_fit;
    intermedss{i} = intermeds;
    energies(i) = E_defs(i) + E_fits(i);
end

[minE, index] = min(energies);
fprintf('Best fit: (%s), with energy: (%f)\n', labels(index), minE);

% tic1 = tic;
% [xs, A, t, E_def, E_fit, intermeds] = fit_model(Ip, ps');
% toc1 = toc(tic1);
% fprintf('Finished fit_model (%.2fs)\n', toc1);

for i = 1 : size(labels,1)
    % Visualize final output
    xs = xss{i};
    A = As{i};
    t = ts{i};
    E_def = E_defs(i);
    E_fit = E_fits(i);
    intermeds = intermedss{i};
    ps = CPs{i};
    
    visualize_model(Ip, xs, ps', A, t, 10);
    suptitle(sprintf('E_tot: %.2f E_def: %.2f E_fit: %.2f', ...
    E_def+E_fit, E_def, E_fit));
    % Animate iterative process
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
    
end
% %% Visualize final output
% visualize_model(Ip, xs, ps', A, t, 10);
% suptitle(sprintf('E_tot: %.2f E_def: %.2f E_fit: %.2f', ...
%     E_def+E_fit, E_def, E_fit));
% %% Animate iterative process
% figure;
% for i=1:length(intermeds)
%     xs_i = intermeds{i}{1};
%     A_i = intermeds{i}{2};
%     t_i = intermeds{i}{3};
%     E_def_i = intermeds{i}{4};
%     E_fit_i = intermeds{i}{5};
%     N_B_i = intermeds{i}{6};
%     visualize_model(Ip, xs_i, ps', A_i, t_i, N_B_i);
%     suptitle(sprintf('Iter %d/%d (E_tot: %.2f E_def: %.2f E_fit: %.2f)', ...
%         i, length(intermeds), ...
%         E_def_i+E_fit_i, E_def_i, E_fit_i));
%     pause
% end
