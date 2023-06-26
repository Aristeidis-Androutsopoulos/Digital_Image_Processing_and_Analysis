% Import the image
f = imread("flower.png");
image = im2double(f);

whos f

min_value = min(image(:));
max_value = max(image(:));

snr_db = 15;

% Convert SNR to linear scale
snr_lin = 10^(snr_db / 10);
image_pow = sum(image(:).^2 / numel(image));
noise_pow = image_pow / snr_lin;

%% Noise Addition and Prep
image_gaussian = imnoise(image, 'gaussian', 0, noise_pow);
% figure, imshow(image_gaussian);
% [r, snr] = psnr(image_gaussian, image);

image_salt_pepper = imnoise(image, 'salt & pepper', 0.25);
% figure, imshow(image_salt_pepper);
% [r, snr] = psnr(image_salt_pepper, image);


moving_average_filter_window_size = 5;
median_filter_window_size = 5;

filter_size_ma = [moving_average_filter_window_size, moving_average_filter_window_size];
filter_size_median = [median_filter_window_size, median_filter_window_size];

%% Gaussian Noise

% Moving Average Filter
moving_average_filter = fspecial('average', filter_size_ma);
image_gaussian_moving_average = imfilter(image_gaussian, moving_average_filter, 'replicate');

% Median Filter
image_gaussian_median = medfilt2(image_gaussian, filter_size_median);

figure, montage({image, image_gaussian, image_gaussian_moving_average, image_gaussian_median})
title('Original Image (Upper Left), With Gaussian Noise (Upper Right), Moving Average Filtered (Lower Left), Median Filtered (Lower Right)');


%% Salt and Pepper Noise

% Moving Average Filter
image_saltpepper_moving_average = imfilter(image_salt_pepper, moving_average_filter, 'replicate');

% Median Filter
image_saltpepper_median = medfilt2(image_salt_pepper, filter_size_median);

figure, montage({image, image_salt_pepper, image_saltpepper_moving_average, image_saltpepper_median})
title('Original Image (Upper Left), With Salt and Pepper Noise (Upper Right), Moving Average Filtered (Lower Left), Median Filtered (Lower Right)');

