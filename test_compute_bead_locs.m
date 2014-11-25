%% Demo usage of compute_bead_locs() and eval_spline()

% ps: defines spline control points
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

nb_k = size(ps,1) - 4 + 1; % nb segments
N_B = nb_k * 10; % Nb. beads.

%% First sample pts on spline
xs = [];
ys = [];

for s=0:0.1:(nb_k-0.1) % avoid endpoint
    xy = eval_spline(s, ps);
    xs = [xs; xy(1)];
    ys = [ys; xy(2)];
end

%% Next compute bead locations
[bs, Ms] = compute_bead_locs(ps, N_B);

%% Plot things
figure;
plot(xs, ys); hold on;
plot(bs(1,:), bs(2,:), 'rx');
legend('Spline', 'Beads');
