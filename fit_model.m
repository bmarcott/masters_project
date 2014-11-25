function [x, T, E_def, E_fit] = fit_model(I, c)
%FIT_MODEL Fits cubic b-spline to image I.
%INPUT
%  matrix I: [h x w]
%  array c: [2 x n]
%    Vector of cubic b-spline control points ('home' locations).
%OUTPUT
%  matrix x:
%  matrix T:
%    Final affine transformation.
nb_c = size(c, 2); % Nb. control points
z = zeros([2*nb_c+6,1]); % Initialize z (spline ctrl pts + affine params)

% I.e. N_0 is relative weight E_def vs E_fit
N_0 = 160; % "Standard" nb. of pixels. [See pg. 31 of thesis]
N_B = sum(sum(I));  % Nb. of inked pixels
pi_n = 0.3; % mixing coef. btwn uniform-noise and gaussian-mixture
EPS = 1e-8;

Bs = []; % Bs(b,:,:) -> 2x2*n matrix mapping ctrl pts to bead locs

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

for iter_var=1:size(anneal_sched, 1)
    B = anneal_sched(iter_var, 1);
    max_inner_iters = anneal_sched(iter_var, 2);
    var_b = anneal_sched(iter_var, 3);
    sigma_b = eye(2) * var_b; % covar. matrix for bead gaussians
    E_tot = E_def() + E_fit();
    for iter_inner=1:max_inner_iters
        %% E-step
        % Evaluate responsiblity of each bead for each pixel
        rs = compute_rs(I, bs, var_b, pi_n);
        %% M-step (part 1)
        % Update control point locations z
        % (1) Compute M_def, b_def
        [M_def, b_def] = compute_Mb_def(sigma_b, A, t);
        % (2) Compute M_fit, b_fit
        [M_fit, b_fit] = compute_Mb_fit(nb_c, B, N_0, N_B, rs, Bs);
        % Solve linear system to update bead locations: Mx = b
        M = (M_def + M_fit);
        b = (b_def + b_fit);
        c_new = (M\b); % stacked (x,y) coords. In image coords.
        %x_new = reshape(x_new, [2, 3])'; % [N x 2]
        %% M-step (part 2)
        % Update affine trans. A,t by minimizing E_def while keeping x_new fixed
        [A, t] = min_E_def(c_new, c_home);

        %% Check stopping criterion
        E_tot_p = E_def() + E_fit();
        delt = E_tot_p - E_tot;
        if abs(delt) <= EPS
            break % Converged
        end
        E_tot = E_tot_p;
    end
end
end

function rs = compute_rs(I, bs, var_b, pi_n)
%INPUT
%  matrix I: [h x w]
%  matrix bs:
%    Stores the bead locations: bs(b,:) -> [x,y]
%  float var_b:
%    Variance of the beads.
%  float pi_n:%
%    Mixing coef. between uniform and gaussian mixture.
%OUTPUT
%  matrix rs:
%    Stores the responsibilities.
%    rs(b,i,j) -> resp. for bead b for pixel at (i,j)
rs = []; % rs(b,i,j) = float responsibility for pixel (i,j)
A = size(I,1)*size(I,2); % Area of image
for b=1:B
    mu_b = bs(b,:)'; % col vec
    for i=1:size(I,1)
        for j=1:size(I,2)
            if (~I(i,j))
                continue
            end
            norm_term = arrayfun(@(mu) (mvnpdf([j,i],mu,var_b*eye(2)) + (pi_n*B)/((1-pi_n)*A)), bs);
            rs(b,i,j) = mvnpdf([j,i], mu_b, var_b*eye(2)) / norm_term;
        end
    end
end
end

function [M_def, b_def] = compute_Mb_def(sigma_0, A, t)
%INPUT
%  matrix sigma_0: [2 x 2]
%    Covariance matrix of the bead gaussians.
%  matrix A: [2 x 2]
%    The affine transformation mapping object frame to image frame.
%  array t: [2 x 1]
%    The translation component of the object->image mapping. [tx, ty].
%OUTPUT
%  matrix M_def: [2 x 2]
%  array b_def: [2 x 1]
sigma_i = inv(A)'*inv(sigma_0)*inv(A);
M_def = inv(sigma_i);
b_def = inv(sigma_i) * (A*h_0 + t);
end

function [M_fit, b_fit] = compute_Mb_fit(nb_c, B, N_0, N_B, rs, Bs)
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
        for ii=1:size(I, 1)
            for jj=1:size(I, 2)
                if ~I(ii,jj)
                    continue
                end
                t2 = t2 + (rs(g,ii,jj)*(Bs(g,1,i)*jj + Bs(g,2,i)*ii));
            end
        end
    end
    b_fit(i) = t1*t2;
end
end
