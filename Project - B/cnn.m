%% Load the MNIST dataset
load ('mnist.mat')

% Training Set 
% Reshaping to [height, width , num_planes, num_samples]

XTrain = reshape(training.images, [28, 28, 1, 60000]);
XTest = reshape(test.images, [28, 28, 1, 10000]);

% Validation Set
YTrain = categorical(training.labels);
YTest = categorical(test.labels);

%% Visualize one image from each class
% Determine the number of unique classes
uniqueClasses = unique(YTrain);

% Create a figure to display one image from each class
figure;

% Loop through each unique class and display one image
for i = 1:numel(uniqueClasses)
    % Find the indices of all images in the class
    classIndices = find(YTrain == uniqueClasses(i));
    
    % Shuffle the indices to get a random order
    shuffledIndices = classIndices(randperm(length(classIndices)));
    
    % Select the first (random) image from the shuffled list
    randomIndex = shuffledIndices(1);
    
    % Extract and display the random image
    subplot(2, 5, i); % Adjust subplot layout as needed
    imshow(XTrain(:, :, 1, randomIndex));
    title(['Class ', num2str(grp2idx(uniqueClasses(i)))]);
end



%% Define your CNN model
layers = [
    % Input Layer
    imageInputLayer([28, 28, 1])
    
    % First Convolution Layer - RELU act. func. and Normalization
    convolution2dLayer(3, 6, 'Padding', 0, 'Stride', 1)
    reluLayer
    batchNormalizationLayer
    
    % Second COnvolution Layer - Avg. Pooling, RELU act. func. and Normalization
    averagePooling2dLayer(2,'Stride',2)
    convolution2dLayer(3, 16, 'Padding', 0, 'Stride', 1)
    reluLayer
    batchNormalizationLayer
    
     % Fully Connected Layer 1st - Avg. Pooling and RELU act. func.
    averagePooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(120)
    reluLayer
    
    % Fully Connected Layer 2nd - RELU act. func.
    fullyConnectedLayer(84)
    reluLayer
    
    % Fully Connected Layer 3rd - SoftMax output act. func.
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer
];

%% Create an options structure
options = trainingOptions('sgdm', ...
    'MaxEpochs', 1, ...
    'MiniBatchSize', 64, ...
    'InitialLearnRate', 0.01, ...
    'ValidationData', {XTest, YTest}, ...
    'ValidationFrequency',30, ...
    'Plots', 'training-progress');

%% Train your CNN model
net = trainNetwork(XTrain, YTrain, layers, options);


% Evaluate the model on the test dataset
YTestPred = classify(net, XTest);
accuracy = sum(YTestPred == YTest) / numel(YTest);
fprintf('Test accuracy: %.2f%%\n', 100 * accuracy);

%% Confusion Matrix
% Convert predicted and true labels to their respective classes
YTestPred = grp2idx(YTestPred); % Convert predicted labels to numeric values
YTestTrue = grp2idx(YTest);     % Convert true labels to numeric values

custom_Confusion_Matrix(YTestPred, YTestTrue)
