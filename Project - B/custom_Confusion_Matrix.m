function custom_Confusion_Matrix(YTestPred, YTestTrue)
    % Determine the number of classes
    numClasses = max(YTestTrue);

    % Initialize the confusion matrix
    confusionMatrix = zeros(numClasses);

    % Populate the confusion matrix
    for i = 1:length(YTestTrue)
        trueClass = YTestTrue(i);
        predClass = YTestPred(i);
        confusionMatrix(trueClass, predClass) = confusionMatrix(trueClass, predClass) + 1;
    end


    % Create a heatmap of the confusion matrix
    figure;
    colormap(parula); % Set the colormap to 'parula'
    % Replace empty cells (NaN or 0) with actual zeros
    confusionMatrix(isnan(confusionMatrix)) = 0;
    heatmap(1:numClasses, 1:numClasses, confusionMatrix, 'ColorLimits', [0, max(confusionMatrix(:))]);

    % Customize the heatmap
    title('Confusion Matrix');
    xlabel('Predicted Class');
    ylabel('True Class');
    colorbar;
end