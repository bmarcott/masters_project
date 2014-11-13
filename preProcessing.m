trainingImagePath = 'C:\Users\Nicholas\Visual Modelling\train-images.idx3-ubyte';
trainingLabelPath = 'C:\Users\Nicholas\Visual Modelling\train-labels.idx1-ubyte';
numToTrain = 20;
offset = 0;

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);
image = imgs(:,:,1);

level = graythresh(image);

subplot(1,2,1)
imshow(image);
title('Before thresholding');

image(image < level) = 1;
image(image ~= 1) = 0;

subplot(1,2,2)
imshow(image);
title('After thresholding');