trainingImagePath = 'C:\Users\Nicholas\Visual Modelling\train-images.idx3-ubyte';
trainingLabelPath = 'C:\Users\Nicholas\Visual Modelling\train-labels.idx1-ubyte';
numToTrain = 20;
offset = 0;

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);
image = imgs(:,:,1);

level = graythresh(image);

subplot(1,3,1)
imshow(image);
title('Before Otsu Thresholding');

image(image < level) = 0;
image(image ~= 0) = 1;

subplot(1,3,2)
imshow(image);
title('After Thresholding');

image = bwmorph(image,'skel',Inf);
subplot(1,3,3)
imshow(image);
title('After Thinning/Skeletonization');