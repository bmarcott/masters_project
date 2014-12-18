function E_fit = compute_E_fit(rs, norm_terms, N_0, N_B);
%COMPUTE_E_FIT Compute E_fit.
%INPUT
%  matrix rs:
%    Stores the responsibilities.
%      rs(b,i,j) -> resp. for bead b for pixel at (i,j)
%  matrix norm_terms:
%    Stores the norm_terms for each pixel:
%      norm_terms(i,j) -> float norm_term
%  int N_0:
%    "standard" nb. of inked pixels.
%  int N_B:
%    Nb. of inked pixels.
%  matrix inked: [N x 2]
%    Stores rows,cols. (i,j).
%    
%OUTPUT
%  float E_fit:
% pg. 31: 
%   E_fit = -(N_0/N_B) * sum(log(f(y_i)) over inked pixels)
% where f(y_i) = (pi_n)/A + (1-pi_n)/B * sum(f_i(y) over bead i)
% and f_i(y) is Norm(mu_i - y, sigma_i)
% TODO: This is wrong! Need the numerators of rs, not denoms
eps = 1e-8;

if 1
    t1 = sum(sum(log(norm_terms(norm_terms > 0.0))));
else
    t1 = 0.0;
    for b=1:size(rs, 1)
        rs_b = squeeze(rs(b,:,:));
        inds = find(norm_terms > eps);
        inds = intersect(inds, inds(find(rs_b(inds) > eps)));
        t1 = t1 + sum(log(rs_b(inds) .* norm_terms(inds)));
    end
end
E_fit = -(N_0/N_B) * t1;
if E_fit < 0
    % TODO: I changed compute_rs() to match what the thesis does (pg.27).
    % But I can't juse tke the denom of rs() anymore, since:
    %   denom(rs()) != f(y)
    % Need to recompute f(y) according to pg. 23
    fprintf('WARNING: E_fit < 0 (%.4f)\n', E_fit);
end

if 0
% Compute it manually
% pg. 23 for f(y)
E_fit2 = 0.0;
inked_inds = find(I ~= 0);
inked_yxs = ind2sub(size(I), inked_inds);
inked_xys(:,1) = inked_yxs(:,2);
inked_xys(:,2) = inked_yxs(:,1);
t = 0;
for b=1:size(rs,1)
    mu_b = bs(:,b);
    t = t + sum(mvnpdf(inked_xys, mu_b', var_b));
end
fy = (pi_n/A) + ((1-pi_n)/size(rs,1)) * t;
end
end
