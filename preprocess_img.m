function Ip = preprocess_img(I)
%PREPROCESS_IMG Summary of this function goes here
%   Detailed explanation goes here
level = graythresh(I);
I(I < level) = 0;
I(I ~= 0) = 1;
Ip = bwmorph(I,'skel',Inf);
end
