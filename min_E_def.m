function [A, t] = min_E_def(c, c_home, LAMBDA)
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
if 0
    %% Buggy way of solving it analytically
    [A, t] = min_v2(c1, c2);
else
    %% Use cvx to solve numerically
    cvx_begin quiet
    %cvx_solver sedumi
    variable A(2, 2)
    variable t(2, 1)
    minimize( compute_E_def(c1, c2, A, t) + LAMBDA*(norm(A) + norm(t)));
    cvx_end
end
end

function [A, t] = min_v2(c, c_home)
%INPUT
%  array c: [2 x N]
%  array c_home: [2 x N]
n = size(c, 2);
% Note: this is wrong! Either the math or code has a bug.
%% Compute M1, b1
M_A = zeros(4,4);
B = c_home(:,1)*c_home(:,1)';
for i=2:n
    B = B + c_home(:,i)*c_home(:,i)';
end
M_A = kron(eye(2, 2), B); % makes blk-diag w/ B's on diag
M_t = n*eye(2);
B1 = zeros(4,4);
for i=1:n
    for j=1:n
        B1 = B1 + kron(inv(M_t)', c_home(:,i)*c_home(:,j)');
    end
end
M1 = -2*pinv(M_A)*B1;
BB1 = zeros(2,2);
for i=1:n
    for j=1:n
        BB1 = BB1 + inv(M_t)*c(:,j)*c_home(:,i)';
    end
    BB1 = BB1 - (c(:,i)*c_home(:,i)');
end
b1 = vec(BB1');

a = M1\b1;
% t* = inv(M_t)*b_t

end
