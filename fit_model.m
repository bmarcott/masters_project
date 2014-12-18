function [xs_est, A, t, E_def, E_fit, intermeds] = fit_model(I, cs_home, anneal_sched, N_0, pi_n, C, lambda_reg, verbose)
%FIT_MODEL Fits cubic b-spline to image I.
%INPUT
%  matrix I: [h x w]
%  array c: [2 x n]
%    Vector of cubic b-spline ctrl pts, 'home' locs. In obj frame.
%OUTPUT
%  matrix xs_est:
%  matrix A, t:
%    Final affine transformation.
%  float E_def, E_fit:
%    The energy terms.
%  cell intermeds: 
%    Outputs the intermediate states:
%      intermeds{i} -> {xs_est, A, t, E_def, E_fit, N_B}
nb_c = size(cs_home, 2); % Nb. control points

% I.e. N_0 is relative weight E_def vs E_fit
%N_0 = 160; % "Standard" nb. of pixels. [See pg. 31 of thesis]
N_I = sum(sum(I));  % Nb. of inked pixels (assumes binary image I)
%pi_n = 0.3; % mixing coef. btwn uniform-noise and gaussian-mixture
EPS = 1e-4;

[inkedRows, inkedCols] = find(I == 1);
inked = [inkedRows, inkedCols]; % [N x 2]

%% Initialize parameters
% Affine trans (A,t) maps object frame to image frame
[A, t] = init_affine(I, cs_home);
% xs_est [2 x N] is ctrl-pt estimates (in image frame)
xs_est = A*cs_home + repmat(t, [1, nb_c]); 
%% Initialize output vals
E_def = nan; E_fit = nan;
E_tots = [];
intermeds = {{xs_est, A, t, nan, nan, 8}};

%% Minimize Energy Function
iter_overall = 0;
for iter_var=1:size(anneal_sched, 1)
    %% Pull out meta-params from annealing schedule
    N_B = anneal_sched(iter_var, 1); % Nb. beads
    max_inner_iters = anneal_sched(iter_var, 2);
    var_b = (anneal_sched(iter_var, 3)) * size(I,1);
    sigma_b = eye(2) * var_b; % covar. matrix for bead gaussians
    [bs, Bs] = compute_bead_locs(xs_est', N_B); % in img frame
    [rs, norm_terms] = compute_rs(I, inked, bs, var_b, pi_n, N_B);
    E_def = compute_E_def(xs_est, cs_home, A, t);
    E_fit = compute_E_fit(rs, norm_terms, pi_n, N_0, N_I);
    E_tot = C*E_def + E_fit;
    E_tots = [E_tots E_tot];
    %% Perform alternating minimization (via EM steps)
    for iter_inner=1:max_inner_iters
                
        [bs, Bs] = compute_bead_locs(xs_est', N_B); % in img frame
        %% E-step
        % Evaluate responsiblity of each bead for each pixel
        [rs, norm_terms] = compute_rs(I, inked, bs, var_b, pi_n, N_B);
        
        E_def_before_1 = compute_E_def(xs_est, cs_home, A, t);
        E_fit_before_1 = compute_E_fit(rs, norm_terms, pi_n, N_0, N_I);
        E_tot_p_before_1 = C*E_def_before_1 + E_fit_before_1;
        
        %% M-step (part 1)
        % Update control point locations xs_est
        % (1) Compute M_def, b_def (see pg. 76)
        [M_def, b_def] = compute_Mb_def(sigma_b, A, t, cs_home);
        % (2) Compute M_fit, b_fit
        [M_fit, b_fit] = compute_Mb_fit(nb_c, N_B, N_0, N_I, rs, Bs, var_b);
        % Solve linear system to update bead locations: Mx = b
        M = (C*M_def + M_fit);
        b = (C*b_def + b_fit);
        cs_new = (M\b); % stacked (x,y) coords. In image coords. M*cs = b
        % Ensure that ctrl points never leave the image bounds
        cs_new(1:2:end) = min(max(cs_new(1:2:end), 0), size(I,2));
        cs_new(2:2:end) = min(max(cs_new(2:2:end), 0), size(I,1));
        xs_est = reshape(cs_new, [2, nb_c]);
        
        [bs, Bs] = compute_bead_locs(xs_est', N_B); % in img frame
        [rs, norm_terms] = compute_rs(I, inked, bs, var_b, pi_n, N_B);
        E_def_after_1 = compute_E_def(xs_est, cs_home, A, t);
        E_fit_after_1 = compute_E_fit(rs, norm_terms, pi_n, N_0, N_I);
        E_tot_p_after_1 = C*E_def_after_1 + E_fit_after_1;
        
        
        %% M-step (part 2)
        % Update affine trans. A,t by minimizing E_def while keeping x_new fixed
        [A, t] = min_E_def(cs_new, cs_home, lambda_reg);
        E_def = compute_E_def(xs_est, cs_home, A, t);
        E_fit = compute_E_fit(rs, norm_terms, pi_n, N_0, N_I);
        %% Update intermeds
        intermeds = [intermeds {{xs_est, A, t, E_def, E_fit, N_B}}];
        %% Check stopping criterion
        E_tot_p = C*E_def + E_fit;
        delt = E_tot_p - E_tot;
        E_tot = E_tot_p;
        E_tots = [E_tots E_tot];
        iter_overall = iter_overall + 1;
        E_compare = [E_def_before_1, E_fit_before_1, E_tot_p_before_1; ...
                     E_def_after_1, E_fit_after_1, E_tot_p_after_1; ...
                     E_def, E_fit, E_tot];

        if verbose
            fprintf('iter=%d [iter_var=%d/%d] iter_inner=%d/%d E_tot: %.2f E_def: %.2f E_fit=%.2f\n', ...
                iter_overall, iter_var, size(anneal_sched, 1), ...
                iter_inner, max_inner_iters, ...
                E_tot, E_def, E_fit);
        end
        if abs(delt) <= EPS
            fprintf('Converged!\n');
            break % Converged
        end
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
sigma_i_inv = inv(AA)'*inv(SS)*inv(AA);
M_def = sigma_i_inv;
h_0 = reshape(cs_home, [size(cs_home,1)*size(cs_home,2), 1]);
b_def = sigma_i_inv * (AA*h_0 + repmat(t, [nb_pts, 1]));
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
t1 = (N_0/((var_b^2)*N_B));        
for i=1:size(M_fit,1)
    for j=1:size(M_fit,2)
        t2 = 0.0;
        for g=1:B
            t3 = sum(sum(rs(g,:,:))); % 41.7% -> 6.4% (remove squeeze)
            t4 = Bs(g,1,i)*Bs(g,1,j) + Bs(g,2,i)*Bs(g,2,j);
            t2 = t2 + (t3*t4);
        end
        M_fit(i,j) = t1*t2;
    end
end
% Compute b_fit (Eq. B.4)
b_fit = zeros([2*nb_c, 1]);
for i=1:size(b_fit,1)
    t2 = 0.0;
    for g=1:B
        if 0
            %% Vectorized, but not any faster than for-loops (oddly)
            foo1 = Bs(g,1,i).*repmat(1:size(rs,3), [size(rs,2),1]);
            foo2 = Bs(g,2,i).*repmat((1:size(rs,2))', [1, size(rs,3)]);
            foo = squeeze(rs(g,:,:)).*(foo1 + foo2);
            t2 = t2 + sum(foo(:));
        else
            %% old way of computing t2
            for ii=1:size(rs, 2)
                for jj=1:size(rs, 3)
                    t2 = t2 + (rs(g,ii,jj)*(Bs(g,1,i)*jj + Bs(g,2,i)*ii)); % 43.2%
                end % 40.5% time spent in this loop
            end
        end
    end
    b_fit(i) = t1*t2;
end
end

function [A, t] = init_affine(I, cs_home)
%INIT_AFFINE Initializes the affine trans. according to pg. 25:
%  Map the bbox around ctrl pts to bbox around img pixels.
%INPUT
%  matrix I: [h x w]
%  matrix cs_home: [2 x N]
%OUTPUT
%  matrix A: [2 x 2]
%  array t: [2 x 1]
%% Compute bounding box around white pixels in I
stats = regionprops(uint8(I), 'BoundingBox');
bbox = stats.BoundingBox; % [x y w h]
bbox(1:2) = bbox(1:2) - 0.5; % x,y are offset by .5
% order: [UL; UR; LL; LR]
bbox_img = [[bbox(1:2)];
    [bbox(1)+bbox(3), bbox(2)];
    [bbox(1), bbox(2) + bbox(4)];
    [bbox(1)+bbox(3), bbox(2)+bbox(4)]];
%% Compute bounding box around spline ctrl points (in obj frame)
x1 = min(cs_home(1,:));
y1 = min(cs_home(2,:));
x2 = max(cs_home(1,:));
y2 = max(cs_home(2,:));
w = x2 - x1; h = y2 - y1;
bbox_obj = [[x1 y1];
    [x1+w, y1];
    [x1, y1+h];
    [x1+w, y1+h]];
tform = cp2tform(bbox_obj, bbox_img, 'affine');
% Note: T is arranged in a weird way (transpose it to get familiar affine)
T = tform.tdata.T'; % transpose
A = T(1:2, 1:2);
t = T(1:2, 3);
end
