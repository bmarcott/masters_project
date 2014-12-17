function E_fit = compute_E_fit(rs, norm_terms, N_0, N_B)
%COMPUTE_E_FIT Compute E_fit.
%INPUT
%  matrix rs:
%  matrix norm_terms:
%  int N_0:
%    "standard" nb. of inked pixels.
%  int N_B:
%    Nb. of inked pixels.
%OUTPUT
%  float E_fit:
% pg. 31: 
%   E_fit = -(N_0/N_B) * sum(log(f(y_i)) over inked pixels)
% where f(y_i) = (pi_n)/A + (1-pi_n)/B * sum(f_i(y) over bead i)
% and f_i(y) is Norm(mu_i - y, sigma_i)
t1 = sum(sum(log(norm_terms(norm_terms > 0.0))));
E_fit = -(N_0/N_B) * t1;
end
