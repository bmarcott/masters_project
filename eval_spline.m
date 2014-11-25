function out = eval_spline(s, ps)
%EVAL_SPLINE Output position s on spline with control points ps.
%INPUT
%  float s:
%    Location param: 0 <= s <= K, where K is the number of 'knots'.
%  matrix ps: [N x 2]
%OUTPUT
%  array out: [2 x 1]
i = floor(s) + 1; % segment number
ps_i = ps(i:i+3, :);
ss = s - floor(s); % effective s for this segment: 0 <= ss <= 1
B_i = [[-1 3 -3 1];
    [3 -6 3 0];
    [-3 0 3 0];
    [1 4 1 0]];
x = (1/6) * [ss^3 ss^2 ss 1] * B_i * ps_i(:,1);
y = (1/6) * [ss^3 ss^2 ss 1] * B_i * ps_i(:,2);
out = [x; y];
end
