function [A, t] = min_E_def(c, c_home)
%MIN_E_DEF Implements the M-step Part 2. Minimize E with respect to the
%affine transformation parameters (A, t).
%INPUT
%  array c: [2*n x 1]
%    Array of spline control points. In image-frame.
%  array c_home: [2*n x 1]
%    Array of 'home' spline control point locations. In object-frame.
%OUTPUT
%  matrix A: [2 x 2]
%    2x2 matrix of affine trans. Maps object frame to image frame.
%  array t: [2 x 1]
%    Translation comp. of affine trans.
n = length(c) / 2;
c1 = reshape(c, [2, n]);
c2 = reshape(c_home, [2, n]);
cvx_begin quiet
variable A(2, 2)
variable t(2, 1)
minimize( compute_E_def(c1, c2, A, t) );
cvx_end
end
