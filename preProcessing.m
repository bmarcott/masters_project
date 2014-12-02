trainingImagePath = '/Users/ansuya/Documents/cs269-proj/mnist/train-images-idx3-ubyte';
trainingLabelPath = '/Users/ansuya/Documents/cs269-proj/mnist/train-labels-idx1-ubyte';
%can set to 60,000 for all training examples and 10,000 for all testing
%examples
numToTrain = 60000;
offset = 0;

[imgs, labels] = readMNIST(trainingImagePath, trainingLabelPath, numToTrain, offset);
for i = 1:numToTrain
image = imgs(:,:,i);

level = graythresh(image);

%subplot(1,3,1)
%imshow(image);
%title('Before Otsu Thresholding');
%mkdir(/Users/ansuya/Documents/cs269-proj/','training-raw');
imwrite(image, strcat('/Users/ansuya/Documents/cs269-proj/training-raw/img', num2str(i), '.jpg'));

image(image < level) = 0;
image(image ~= 0) = 1;

%subplot(1,3,2)
%imshow(image);
%title('After Thresholding');

image = bwmorph(image,'skel',Inf);
%subplot(1,3,3)
%imshow(image);
%title('After Thinning/Skeletonization');
%mkdir(/Users/ansuya/Documents/cs269-proj/','training-processed');
imwrite(image, strcat('/Users/ansuya/Documents/cs269-proj/training-processed/img', num2str(i), '.jpg'));
end