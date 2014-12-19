

function [bs, Ms] = compute_bead_locs(cs, N_B)
%COMPUTE_BEAD_LOCS Compute the bead locations along spline cs.
%INPUT
%  array cs: [N x 2]
%  int N_B:
%    Nb. of beads.
%OUTPUT
%  matrix bs: [2 x N_B]
%    Contains the bead locations. Cols are xy coords.
%  matrix Ms: [N_B x 2 x 2n]
%    Contains the [2 x 2n] matrices that map spline ctrl pts to bead locs.
n = size(cs, 1);
bs = [];
Ms = zeros(N_B, 2, 2*n);


nb_k = size(cs,1) - 4 + 1; % nb segments

xs = [];
ys = [];
parms = [];

totalArcLength = 0;
for s=0:0.001:(nb_k-0.1) % avoid endpoint
    xy = eval_spline(s, cs);
    xs = [xs; xy(1)];
    ys = [ys; xy(2)];
    parms = [parms;s];
    if s > 0
        x1 = xs(end-1);
        x2 = xs(end);
        y1 = ys(end-1);
        y2 = ys(end);
        totalArcLength = totalArcLength + sqrt((x2-x1)^2+(y2-y1)^2);
    end
end

incArcLength = totalArcLength/(N_B+1);

curArcLength = 0;
beadParms = [];
for i = 2 : size(xs,1);
    
    x1 = xs(i-1);
    x2 = xs(i);
    y1 = ys(i-1);
    y2 = ys(i);
    curArcLength = curArcLength + sqrt((x2-x1)^2+(y2-y1)^2);
    if curArcLength >= incArcLength
        beadParms = [beadParms parms(i)];
        curArcLength = 0;
    end
    
    
end


for b=1:N_B
    [xy, M] = compute_bead_loc2(beadParms(b), cs, N_B);
    bs = [bs xy];
    Ms(b, :, :) = M;
end
end

function [out, M] = compute_bead_loc2(s, ps, N_B)

N = size(ps, 1); % nb. ctrl points

seg_cur = floor(s+1);
ss = s+1 - seg_cur; % effective s for this segment

B_i = [[-1 3 -3 1];
    [3 -6 3 0];
    [-3 0 3 0];
    [1 4 1 0]];

S = [[[ss^3 ss^2 ss 1] zeros(1,4)];
    [zeros(1,4) [ss^3 ss^2 ss 1]]];
B = [[B_i zeros(4,4)];
    [zeros(4,4) B_i]];
% G*p extracts the 4 relevant ctrl pts, in [x1;y1;...] stacks
G = zeros(8, 2*N);
for i=1:size(G,1)
    for j=1:size(G,2)
        if (j == (i + 2*seg_cur - 2))
            G(i,j) = 1;
        end
    end
end
% PERM*(Gp) reorders the [x1;y1;...] stacks to: [x1;x2;...;y1;y2;...]
PERM = zeros(8,8);
for i=1:4
    PERM(i, 2*(i-1)+1) = 1;
end
for i=5:8
    PERM(i, 2*(i-4-1)+2) = 1;
end
M = (1/6)*S*B*PERM*G;
p = reshape(ps', [2*N, 1]); % stack pts [x1;y1;x2;y2;...]
out = M*p;

end

% function [bs, Ms] = compute_bead_locs(cs, N_B)
% %COMPUTE_BEAD_LOCS Compute the bead locations along spline cs.
% %INPUT
% %  array cs: [N x 2]
% %  int N_B:
% %    Nb. of beads.
% %OUTPUT
% %  matrix bs: [2 x N_B]
% %    Contains the bead locations. Cols are xy coords.
% %  matrix Ms: [N_B x 2 x 2n]
% %    Contains the [2 x 2n] matrices that map spline ctrl pts to bead locs.
% n = size(cs, 1);
% bs = [];
% Ms = zeros(N_B, 2, 2*n);
% for b=1:N_B
%     [xy, M] = compute_bead_loc2(b, cs, N_B);
%     bs = [bs xy];
%     Ms(b, :, :) = M;
% end
% end
% 
% function [out, M] = compute_bead_loc(b, ps, N_B)
% COMPUTE_BEAD_LOC Compute location of bead b on spline ps. Express the
% operation as a matrix-vector operation on *all* ctrl points.
% Assumes N_B is divisible by nb segments
% Also outputs the [2 x 2n] matrix M that maps ctrl pts to bead locs.
% INPUT
%  int b:
%    Bead number.
%  matrix ps: [N x 2]
%    The spline control points.
%  int N_B:
%    Nb. of beads total. Must be divisible by nb. segments.
% OUTPUT
%  array xy: [2 x 1]
%    The location for bead b on the spline.
%  matrix M: [2 x 2n]
%    The matrix that maps spline ctrl pts to the b-th bead location.
% N = size(ps, 1); % nb. ctrl points
% nb_k = size(ps,1) - 4 + 1; % nb segments
% beads_per_seg = N_B / nb_k;
% seg_cur = ceil(b / beads_per_seg);
% ss = mod(b, beads_per_seg) / beads_per_seg; % effective s for this segment
% 
% B_i = [[-1 3 -3 1];
%     [3 -6 3 0];
%     [-3 0 3 0];
%     [1 4 1 0]];
% 
% S = [[[ss^3 ss^2 ss 1] zeros(1,4)];
%     [zeros(1,4) [ss^3 ss^2 ss 1]]];
% B = [[B_i zeros(4,4)];
%     [zeros(4,4) B_i]];
% G*p extracts the 4 relevant ctrl pts, in [x1;y1;...] stacks
% G = zeros(8, 2*N);
% for i=1:size(G,1)
%     for j=1:size(G,2)
%         if (j == (i + 2*seg_cur - 2))
%             G(i,j) = 1;
%         end
%     end
% end
% PERM*(Gp) reorders the [x1;y1;...] stacks to: [x1;x2;...;y1;y2;...]
% PERM = zeros(8,8);
% for i=1:4
%     PERM(i, 2*(i-1)+1) = 1;
% end
% for i=5:8
%     PERM(i, 2*(i-4-1)+2) = 1;
% end
% M = (1/6)*S*B*PERM*G;
% p = reshape(ps', [2*N, 1]); % stack pts [x1;y1;x2;y2;...]
% out = M*p;
% end
