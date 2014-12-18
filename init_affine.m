function [A, t] = init_affine(I, cs_home)
%INIT_AFFINE Initializes the affine trans. according to pg. 25:
%  Map the bbox around ctrl pts to bbox around img pixels.
%INPUT
%  matrix I: [h x w]
%  matrix cs_home: [2 x N]
%OUTPUT
%  matrix A: [2 x 2]
%  array t: [2 x 1]
%% Compute bounding box around white pixels in I
stats = regionprops(uint8(I), 'BoundingBox');
bbox = stats.BoundingBox; % [x y w h]
bbox(1:2) = bbox(1:2) - 0.5; % x,y are offset by .5
% order: [UL; UR; LL; LR]
bbox_img = [[bbox(1:2)];
    [bbox(1)+bbox(3), bbox(2)];
    [bbox(1), bbox(2) + bbox(4)];
    [bbox(1)+bbox(3), bbox(2)+bbox(4)]];
%% Compute bounding box around spline ctrl points (in obj frame)
x1 = min(cs_home(1,:));
y1 = min(cs_home(2,:));
x2 = max(cs_home(1,:));
y2 = max(cs_home(2,:));
w = x2 - x1; h = y2 - y1;
bbox_obj = [[x1 y1];
    [x1+w, y1];
    [x1, y1+h];
    [x1+w, y1+h]];
tform = cp2tform(bbox_obj, bbox_img, 'affine');
% Note: T is arranged in a weird way (transpose it to get familiar affine)
T = tform.tdata.T'; % transpose
A = T(1:2, 1:2);
t = T(1:2, 3);
end
