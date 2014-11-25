function [xs_est, A, t, E_def, E_fit] = fit_model(I, cs_home)
%FIT_MODEL Fits cubic b-spline to image I.
%INPUT
%  matrix I: [h x w]
%  array c: [2 x n]
%    Vector of cubic b-spline control points ('home' locations).
%OUTPUT
%  matrix xs_est:
%  matrix A, t:
%    Final affine transformation.
%  float E_def, E_fit:
%    The energy terms.
nb_c = size(cs_home, 2); % Nb. control points

% I.e. N_0 is relative weight E_def vs E_fit
N_0 = 160; % "Standard" nb. of pixels. [See pg. 31 of thesis]
N_I = sum(sum(I));  % Nb. of inked pixels (assumes binary image I)
pi_n = 0.3; % mixing coef. btwn uniform-noise and gaussian-mixture
EPS = 1e-8;

anneal_sched = [...
    %[int NB_BEADS, int MAX_INNER_ITERS, float VAR]
    [8, 5, 0.04]; ...
    [15, 3, 0.02]; ...
    [15, 3, 0.01]; ...
    [15, 35, 0.008]; ...
    [30, 35, 0.0025]; ...
    [60, 45, 0.0006]; ...
    [60, 2, 0.0006]; ...
    ];

MAX_ITERS = 3;

%% Initialize parameters
% Affine trans (A,t) maps object frame to image frame
A = eye(2); % TODO: Actually init. this trans.
t = [0; 0];
% xs_est [2 x N] is ctrl-pt estimates (in image frame)
xs_est = A*cs_home + repmat(t, [1, nb_c]); 

%% Initialize output vals
E_def = nan; E_fit = nan;
E_tots = [];
for iter_var=1:size(anneal_sched, 1)
    %% Pull out meta-params from annealing schedule
    N_B = anneal_sched(iter_var, 1); % Nb. beads
    max_inner_iters = anneal_sched(iter_var, 2);
    var_b = anneal_sched(iter_var, 3);
    sigma_b = eye(2) * var_b; % covar. matrix for bead gaussians
    [bs, Bs] = compute_bead_locs(xs_est', N_B); % in img frame
    [rs, norm_terms] = compute_rs(I, bs, var_b, pi_n, N_B);
    E_def = compute_E_def(xs_est, cs_home, A, t);
    E_fit = compute_E_fit(rs, norm_terms, N_0, N_B);
    E_tot = E_def + E_fit;
    E_tots = [E_tots E_tot];
    %% Perform alternating minimization (via EM steps)
    for iter_inner=1:max_inner_iters
        if iter_inner > MAX_ITERS
            break % DEBUG
        end
        fprintf('[iter_var=%d/%d] iter_inner=%d/%d\n', ...
            iter_var, size(anneal_sched, 1), ...
            iter_inner, max_inner_iters);
        fprintf('E_tots: [%s]\n', sprintf('%.4f ', E_tots));
        [bs, Bs] = compute_bead_locs(xs_est', N_B); % in img frame
        %% E-step
        % Evaluate responsiblity of each bead for each pixel
        [rs, norm_terms] = compute_rs(I, bs, var_b, pi_n, N_B);
        %% M-step (part 1)
        % Update control point locations xs_est
        % (1) Compute M_def, b_def
        [M_def, b_def] = compute_Mb_def(sigma_b, A, t, cs_home);
        % (2) Compute M_fit, b_fit
        [M_fit, b_fit] = compute_Mb_fit(nb_c, N_B, N_0, N_I, rs, Bs, var_b);
        % Solve linear system to update bead locations: Mx = b
        M = (M_def + M_fit);
        b = (b_def + b_fit);
        cs_new = (M\b); % stacked (x,y) coords. In image coords.
        xs_est = reshape(cs_new, [2, nb_c]);
        %% M-step (part 2)
        % Update affine trans. A,t by minimizing E_def while keeping x_new fixed
        [A, t] = min_E_def(cs_new, cs_home);

        %% Check stopping criterion
        E_def = compute_E_def(xs_est, cs_home, A, t);
        E_fit = compute_E_fit(rs, norm_terms, N_0, N_B);
        E_tot_p = E_def + E_fit;
        delt = E_tot_p - E_tot;
        if abs(delt) <= EPS
            break % Converged
        end
        E_tot = E_tot_p;
        E_tots = [E_tots E_tot];
    end
end
end

function [M_def, b_def] = compute_Mb_def(sigma_0, A, t, cs_home)
%INPUT
%  matrix sigma_0: [2 x 2]
%    Covariance matrix of the bead gaussians.
%  matrix A: [2 x 2]
%    The affine transformation mapping object frame to image frame.
%  array t: [2 x 1]
%    The translation component of the object->image mapping. [tx, ty].
%  matrix cs_home: [2 x N]
%    Home locations of spline ctrl pts. In obj frame.
%OUTPUT
%  matrix M_def: [2 x 2]
%  array b_def: [2 x 1]
nb_pts = size(cs_home, 2);
AA = kron(eye(nb_pts, nb_pts), A); % makes blk-diag w/ A's on diag
SS = kron(eye(nb_pts, nb_pts), sigma_0);
sigma_i = inv(AA)'*inv(SS)*inv(AA);
M_def = inv(sigma_i);
h_0 = reshape(cs_home, [size(cs_home,1)*size(cs_home,2), 1]);
b_def = inv(sigma_i) * (AA*h_0 + repmat(t, [nb_pts, 1]));
end

function [M_fit, b_fit] = compute_Mb_fit(nb_c, B, N_0, N_B, rs, Bs, var_b)
%INPUT
%  int nb_c:
%    Nb. of control points.
%  int B:
%    Nb. of beads.
%  float N_0, N_B:
%    (N_0/N_B) determines relative importance of inked pixels.
%  matrix rs:
%    Contains the responsibilities:
%      rs(bead,i,j) -> Resp. of bead at location (i,j).
%  matrix Bs:
%    Contains the Bg [2 x 2*n] matrices for each bead:
%      Bs(bead,:,:) -> [2 x 2n] matrix for bead
%  
% Compute M_fit (Eq. B.3)
M_fit = zeros([2*nb_c, 2*nb_c]);
t1 = (N_0/((var_b^2)*B));        
for i=1:size(M_fit,1)
    for j=1:size(M_fit,2)
        t2 = 0.0;
        for g=1:B
            t3 = sum(sum(squeeze(rs(g,:,:))));
            t4 = Bs(g,1,i)*Bs(g,1,j) + Bs(g,2,i)*Bs(g,2,j);
            t2 = t2 + (t3*t4);
        end
        M_fit(i,j) = t1*t2;
    end
end
% Compute b_fit (Eq. B.4)
b_fit = zeros([2*nb_c, 1]);
t1 = (N_0/((var_b^2)*B));
for i=1:size(b_fit,1)
    t2 = 0.0;
    for g=1:B
        for ii=1:size(rs, 2)
            for jj=1:size(rs, 3)
                t2 = t2 + (rs(g,ii,jj)*(Bs(g,1,i)*jj + Bs(g,2,i)*ii));
            end
        end
    end
    b_fit(i) = t1*t2;
end
end
