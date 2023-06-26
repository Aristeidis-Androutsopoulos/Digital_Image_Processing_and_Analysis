% Import the image
road_1 = imread("dark_road_1.jpg");
road_2 = imread("dark_road_2.jpg");
road_3 = imread("dark_road_3.jpg");

% whos road_1
% whos road_2
% whos road_3

%% Original Images
road_1 = im2double(road_1);
road_2 = im2double(road_2);
road_3 = im2double(road_3);

% figure, montage({road_1, road_2, road_3})

% Histograms of Original Images
road_1_hist = imhist(road_1);
road_2_hist = imhist(road_2);
road_3_hist = imhist(road_3);

%% Global Histogram Equalized Images
road_1_eq = histeq(road_1);
road_2_eq = histeq(road_2);
road_3_eq = histeq(road_3);

% figure, montage({road_1_eq, road_2, road_3})

road_1_eq_hist = imhist(road_1_eq);
road_2_eq_hist = imhist(road_2_eq);
road_3_eq_hist = imhist(road_3_eq);

%% Local Histogram Equalized Images
size = 64;
window_size = [size, size];

local_histogram_fun = @(block) histeq(block.data);
road_1_local_eq = blockproc(road_1, window_size, local_histogram_fun);
road_1_local_eq = rescale(road_1_local_eq);

road_2_local_eq = blockproc(road_2, window_size, local_histogram_fun);
road_2_local_eq = rescale(road_2_local_eq);

road_3_local_eq = blockproc(road_3, window_size, local_histogram_fun);
road_3_local_eq = rescale(road_3_local_eq);

road_1_local_eq_hist = imhist(road_1_local_eq);
road_2_local_eq_hist = imhist(road_2_local_eq);
road_3_local_eq_hist = imhist(road_3_local_eq);


%% Plot all different versions of the images
show_hist(road_1, road_1_hist, road_1_eq, road_1_eq_hist, road_1_local_eq, road_1_local_eq_hist);
show_hist(road_2, road_2_hist, road_2_eq, road_2_eq_hist, road_2_local_eq, road_2_local_eq_hist);
show_hist(road_3, road_3_hist, road_3_eq, road_3_eq_hist, road_3_local_eq, road_3_local_eq_hist);

function show_hist(image, hist, image2, hist2, image3, hist3)
    figure;
    subplot(3, 2, 1);
    imshow(image);
    title('Original Image');
    
    subplot(3, 2, 2);
    bar(hist);
    title('Image Histogram');
    xlabel('Pixel Intensity');
    ylabel('Frequency');
    
    subplot(3, 2, 3);
    imshow(image2);
    title('Global EQ Image');
    
    subplot(3, 2, 4);
    bar(hist2);
    title('Image Histogram');
    xlabel('Pixel Intensity');
    ylabel('Frequency');
    
    subplot(3, 2, 5);
    imshow(image3);
    title('Local EQ Image');
    
    subplot(3, 2, 6);
    bar(hist3);
    title('Image Histogram');
    xlabel('Pixel Intensity');
    ylabel('Frequency');

    % Adjust the subplot arrangement and spacing
    set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
end