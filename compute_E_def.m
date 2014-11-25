function E_def = compute_E_def(cs, cs_home, A, t)
%COMPUTE_E_DEF Compute the E_deform term.
%INPUT
%  matrix cs: [2 x N]
%    Estimate for spline ctrl pts. In image frame.
%  matrix cs_home: [2 x N]
%    Home locations. In object frame.
%  matrix A: [2 x 2]
%    Maps object frame to image frame.
%  array t: [2 x 1]
%OUTPUT
%  float E_def
n = size(cs, 2);
E_def = sum(sum((cs - (A*cs_home + repmat(t,[1,n]))).^2, 1));
end
