%% The function returns the correct size of the feature Vector but has different values than the built in function and thus it does not work properly
% Mistake seems to be part of the binning process.
function features = customHOG(img, cellSize, blockSize, numBins)
    % Calculate gradients using Sobel filters
    [dx, dy] = imgradientxy(image, 'Sobel');

    % Calculate gradient magnitude and orientation
    magnitude = sqrt(dx.^2 + dy.^2);
    orientation = atan2(dy, dx);

    % Determine the size of the HOG grid
    [imageHeight, imageWidth] = size(image);

    % Initialize the HOG feature vector
    features = [];
    
    block = blockSize(1) * cellSize(1);
    jump = size(image,1) - block;
    
    
    
    %% Iterate through blocks and compute HOG features with block interleaving
    for y = 1:jump:2*jump
        for x = 1:jump:2*jump
            % Initialize block histogram
            blockHistogram = [];

            % Iterate through cells within the block
            for by = 0:1
                for bx = 0:1
                    % Define cell boundaries
                    xStart = x + bx * cellSize(1);
                    xEnd = xStart + cellSize(1) - 1;
                    yStart = y + by * cellSize(2);
                    yEnd = yStart + cellSize(2) - 1;
                    
                    % Select gradients and orientations within the cell
                    cellMagnitude = magnitude(yStart:yEnd, xStart:xEnd);
                    cellOrientation = orientation(yStart:yEnd, xStart:xEnd);
                    
                    histr=zeros(1,numBins);
                    
                    for k = 1:8
                        for l = 1:8
                            
                            alpha = cellOrientation(k,l);
                        
                            % Binning Process (Bi-Linear Interpolation)
                            if alpha>10 && alpha<=30
                                histr(1)=histr(1)+ cellMagnitude(k,l)*(30-alpha)/20;
                                histr(2)=histr(2)+ cellMagnitude(k,l)*(alpha-10)/20;
                            elseif alpha>30 && alpha<=50
                                histr(2)=histr(2)+ cellMagnitude(k,l)*(50-alpha)/20;                 
                                histr(3)=histr(3)+ cellMagnitude(k,l)*(alpha-30)/20;
                            elseif alpha>50 && alpha<=70
                                histr(3)=histr(3)+ cellMagnitude(k,l)*(70-alpha)/20;
                                histr(4)=histr(4)+ cellMagnitude(k,l)*(alpha-50)/20;
                            elseif alpha>70 && alpha<=90
                                histr(4)=histr(4)+ cellMagnitude(k,l)*(90-alpha)/20;
                                histr(5)=histr(5)+ cellMagnitude(k,l)*(alpha-70)/20;
                            elseif alpha>90 && alpha<=110
                                histr(5)=histr(5)+ cellMagnitude(k,l)*(110-alpha)/20;
                                histr(6)=histr(6)+ cellMagnitude(k,l)*(alpha-90)/20;
                            elseif alpha>110 && alpha<=130
                                histr(6)=histr(6)+ cellMagnitude(k,l)*(130-alpha)/20;
                                histr(7)=histr(7)+ cellMagnitude(k,l)*(alpha-110)/20;
                            elseif alpha>130 && alpha<=150
                                histr(7)=histr(7)+ cellMagnitude(k,l)*(150-alpha)/20;
                                histr(8)=histr(8)+ cellMagnitude(k,l)*(alpha-130)/20;
                            elseif alpha>150 && alpha<=170
                                histr(8)=histr(8)+ cellMagnitude(k,l)*(170-alpha)/20;
                                histr(9)=histr(9)+ cellMagnitude(k,l)*(alpha-150)/20;
                            elseif alpha>=0 && alpha<=10
                                histr(1)=histr(1)+ cellMagnitude(k,l)*(alpha+10)/20;
                                histr(9)=histr(9)+ cellMagnitude(k,l)*(10-alpha)/20;
                            elseif alpha>170 && alpha<=180
                                histr(9)=histr(9)+ cellMagnitude(k,l)*(190-alpha)/20;
                                histr(1)=histr(1)+ cellMagnitude(k,l)*(alpha-170)/20;
                            end
                        end
                    end
                    
                    % Concatenate the block histogram to the feature vector
                    blockHistogram = [blockHistogram histr];
                end
            end
            
            % L2 normalization within each block
            blockNorm = norm(blockHistogram);
            blockFeatures = blockHistogram / blockNorm;
            features = [features blockFeatures];
        end
    end
end

%% Function to create the cell Histogram
function cellHistogram = computeCellHistogram(cellMagnitude, cellOrientation, numBins)
    cellHistogram = zeros(1, numBins);
    for bin = 1:numBins
        binRange = [(bin - 1) * (180 / numBins), bin * (180 / numBins)];
        binMask = (cellOrientation >= binRange(1)) & (cellOrientation < binRange(2));
        cellHistogram(bin) = sum(cellMagnitude(binMask));
    end
end