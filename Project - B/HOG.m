%% Load the MNIST dataset
load ('mnist.mat')

% Training Set 
% Reshaping to [height, width , num_planes, num_samples]

XTrain = reshape(training.images, [28, 28, 1, 60000]);
XTest = reshape(test.images, [28, 28, 1, 10000]);

% Validation Set
YTrain = categorical(training.labels);
YTest = categorical(test.labels);


%% Visualize the Different Windows

cellSize = [8 8];
blockSize = [2 2];

img = XTrain(:, :, 1, 115);

% Extract HOG features and HOG visualization
[hog_2x2, vis2x2] = extractHOGFeatures(img,'CellSize',[2 2], 'BlockSize', blockSize);
[hog_4x4, vis4x4] = extractHOGFeatures(img,'CellSize',[4 4], 'BlockSize', blockSize);
[hog_8x8, vis8x8] = extractHOGFeatures(img,'CellSize',[8 8], 'BlockSize', blockSize);
[hog_14x14, vis14x14] = extractHOGFeatures(img,'CellSize',[14 14], 'BlockSize', blockSize);

% % Show the original image
% figure; 
% subplot(2,4,1:4); imshow(img);
% 
% 
% % Visualize the HOG features
% subplot(2,4,5); 
% imshow(img); 
% hold on;
% plot(vis2x2, 'Color', 'green'); 
% title({'CellSize = [2 2]'; ['Length = ' num2str(length(hog_2x2))]});
% % Visualize the HOG features
% subplot(2,4,6);
% imshow(img); 
% hold on;
% plot(vis4x4, 'Color', 'green'); 
% title({'CellSize = [4 4]'; ['Length = ' num2str(length(hog_4x4))]});
% % Visualize the HOG features
% subplot(2,4,7);  
% imshow(img); 
% hold on;
% plot(vis8x8, 'Color', 'green'); 
% title({'CellSize = [8 8]'; ['Length = ' num2str(length(hog_8x8))]});
% % Visualize the HOG features
% subplot(2,4,8);  
% imshow(img); 
% hold on;
% plot(vis14x14, 'Color', 'green'); 
% title({'CellSize = [16 16]'; ['Length = ' num2str(length(hog_16x16))]});

%% Create the HOG Vectors for each image

% Extract HOG features for training data
hogFeatureSize = length(hog_8x8);
numImagesTrain = size(XTrain, 4);
numImagesTest = size(XTest, 4);

hog_training_features = zeros(numImagesTrain, hogFeatureSize, 'single'); % Adjust the feature dimension as needed
custom_hog_training_features = zeros(numImagesTrain, hogFeatureSize, 'single');
hog_test_features = zeros(numImagesTest, hogFeatureSize, 'single'); % Adjust the feature dimension as needed
custom_hog_test_features = zeros(numImagesTest, hogFeatureSize, 'single');



for i = 1:numImagesTrain
    if (i<=numImagesTest)
        img_test = XTest(:, :, 1, i);
        img_test = imbinarize(img_test);
        hog_test_features(i, :) = extractHOGFeatures(img_test, 'CellSize', cellSize, 'BlockSize', blockSize);
%         custom_hog_test_features(i, :) = customHOG(img_test, cellSize, blockSize, 9);
    end
    
    img_train = XTrain(:, :, 1, i);    
    img_train = imbinarize(img_train);
    hog_training_features(i, :) = extractHOGFeatures(img_train, 'CellSize', cellSize, 'BlockSize', blockSize);
%     custom_hog_training_features(i, :) = customHOG(img_train, cellSize, blockSize, 9);
    i

end

disp("done")

%% Training and Prediction

% Train an SVM classifier
svmClassifier = fitcecoc(hog_training_features, YTrain);
% custom_svmClassifier = fitcecoc(custom_hog_training_features, YTrain);

% Predict using the SVM classifier
[predictedLabels,NegLoss,PBScore] = predict(svmClassifier, hog_test_features);
% custom_predictedLabels = predict(custom_svmClassifier, hog_test_features);

% Evaluate the SVM classifier
accuracy = sum(predictedLabels == YTest) / numel(YTest);
fprintf('Test accuracy: %.2f%%\n', 100 * accuracy);

% accuracy = sum(custom_predictedLabels == YTest) / numel(YTest);
% fprintf('Test accuracy: %.2f%%\n', 100 * accuracy);

%% Some Prints and Visualizations
% idx = randsample(size(YTest,1),10,1);
% 
% table(YTest(idx),predictedLabels(idx),...
%     'VariableNames',{'TrueLabel','PredLabel'})
% 
% NegLoss(idx,:)


%% Confusion Matrix
YTestPred = grp2idx(predictedLabels); % Convert predicted labels to numeric values
YTestTrue = grp2idx(YTest);     % Convert true labels to numeric values

custom_Confusion_Matrix(YTestPred, YTestTrue)
