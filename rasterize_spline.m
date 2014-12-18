function Ishow = rasterize_spline(xs, size_out)
%RASTERIZE_SPLINE Summary of this function goes here
%INPUT
%  matrix xs: [2 x N]
%    The spline control points in the image frame.
%  array size_out: [h x w]
%    The desired size of the output image.
%OUTPUT
%  
nb_k = size(xs, 2) - 4 + 1; % nb segments
Ishow = zeros(size_out);
%% Draw xs spline onto Ishow
xys = []; % [2 x N]
for s=0:0.1:(nb_k-0.01)
    xy = eval_spline(s, xs');
    xy = round(xy);
    xy(1) = min(max(1, xy(1)), size(Ishow, 1));
    xy(2) = min(max(1, xy(2)), size(Ishow, 2));
    Ishow(xy(2),xy(1)) = 1.0;
end
end
