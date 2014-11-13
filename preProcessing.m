
imageDir = 'C:\Users\Nicholas\Machine Learning\project 2\face_data\newface16\';
imageFiles = dir([imageDir '*.bmp']);     

image = imread([imageDir, imageFiles(1).name]);
level = graythresh(image)*255;

imshow(image);
image(image < level) = 0;
image(image ~= 0) = 255;

figure;
imshow(image);