function [x, T, E_def, E_fit] = fit_model(I, c)
%FIT_MODEL Fits cubic b-spline to image I.
%INPUT
%  matrix I: [h x w]
%  array c:
%    Vector of cubic b-spline control points ('home' locations).
%OUTPUT
%  matrix x:
%  matrix T:
%    Final affine transformation.
nb_c = size(c, 1); % Nb. control points
z = zeros([2*nb_c+6,1]); % Initialize z (spline ctrl pts + affine params)

% I.e. N_0 is relative weight E_def vs E_fit
N_0 = 160; % "Standard" nb. of pixels. [See pg. 31 of thesis]
pi_n = 0.3; % mixing coef. btwn uniform-noise and gaussian-mixture
nb_var_steps = 30; % nb. variance annealing steps
EPS = 1e-8;

Bs = []; % Bs(b,:,:) -> 2x2*n matrix mapping ctrl pts to bead locs

for iter_var=1:nb_var_steps
    % TODO: Set B, var_b according to schedule
    B = 8;
    var_b = 4;
    E_tot = E_def() + E_fit();
    for iter_inner=1:max_iters
        %% E-step
        % Evaluate responsiblity of each bead for each pixel
        rs = compute_rs(I, bs, var_b, pi_n);
        %% M-step (part 1)
        % Update control point locations z
        % (1) Compute M_def, b_def
        % Compute M_def
        
        % (2) Compute M_fit, b_fit
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
        %% M-step (part 2)
        % Update affine transformation T
        
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
