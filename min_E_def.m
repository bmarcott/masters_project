function [A, t] = min_E_def(c, c_home)
%INPUT
%  array c: [2*n x 1]
%    Array of spline control points. In image-frame.
%  array c_home: [2*n x 1]
%    Array of 'home' spline control point locations. In object-frame.
%OUTPUT
%  matrix A: [2 x 2]
%    2x2 matrix of affine trans. Maps image frame to object frame.
%  array t: [2 x 1]
%    Translation comp. of affine trans. Maps img frame -> obj frame.
n = length(c) / 2;
c1 = reshape(c, [2, n]);
c2 = reshape(c_home, [2, n]);
cvx_begin
variable A(2, 2)
variable t(2, 1)
minimize( sum(sum((A*c1 - c2).^2, 1)) )
cvx_end
end
