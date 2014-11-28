function visualize_model(I, xs, ps, A, t, N_B)
%VISUALIZE_MODEL Visualize the model fit.
%INPUT
%  matrix I: [h x w]
%    Assume I is binary
%  matrix xs: [2 x N]
%    The spline control points in img frame.
%  matrix ps: [2 x N]
%    The spline 'home' ctrl pts in obj frame.
%  matrix A: [2 x 2]
%  array t: [2 x 1]
%    The affine transformation mapping obj -> img frame.
%  int N_B:
%    Number of beads.
nb_k = size(xs, 2) - 4 + 1; % nb segments
Ishow = zeros([size(I,1), size(I,2), 3]);
Ishow(:,:,1) = I;
Ishow(:,:,2) = I;
Ishow(:,:,3) = I;
Ishow2 = Ishow;
%% Draw xs spline onto Ishow, Ishow2
xys = []; % [2 x N]
for s=0:0.1:(nb_k-0.01)
    xy = eval_spline(s, xs');
    xy = round(xy);
    xys = [xys xy];
    Ishow(xy(2),xy(1),:) = [1.0; 0; 0];
    Ishow2(xy(2),xy(1),:) = [1.0; 0; 0];
end
%% Draw ps spline onto Ishow2
xys2 = []; % [2 x N]
for s=0:0.1:(nb_k-0.01)
    xy = eval_spline(s, ps');
    xy = round(A*xy + t);
    xys2 = [xys2 xy];
    Ishow2(xy(2),xy(1),:) = [0; 0; 1.0];
    Ishow2(xy(2),xy(1),:) = [0; 0; 1.0];
end

%% Show canonical frame
Ishow_canon = zeros([size(I,1), size(I,2), 3]);
xs_canon = inv(A)*xs - repmat(t, [1, size(xs, 2)]);
for s=0:0.1:(nb_k-0.01)
    xy = eval_spline(s, xs_canon');
    xy = round(xy);
    xy(1) = min(max(1, xy(1)), size(I, 1));
    xy(2) = min(max(1, xy(2)), size(I, 2));
    Ishow_canon(xy(2), xy(1), :) = [1.0; 0; 0];
    xy = eval_spline(s, ps');
    xy = round(xy);
    Ishow_canon(xy(2), xy(1), :) = [0; 0; 1.0];
end
bs = compute_bead_locs(xs', N_B);
subplot(1,4,1);
imshow(I);
title('I');
subplot(1,4,2);
imshow(Ishow); hold on;
plot(bs(1,:), bs(2,:), 'g+');
title('Spline Est (in red)');
subplot(1,4,3);
imshow(Ishow2);
title('Home spline (blue)');
subplot(1,4,4);
imshow(Ishow_canon);
title('Canonical frame');
end
