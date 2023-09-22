% Import the image
f = imread("factory.jpg");
f = rgb2gray(f);
factory = im2double(f);

whos f

% 2D delta impulse signal
impulse = zeros(31, 31);
impulse(16, 16) = 1;

%% Gauss Smoothing Kernel Filter

std = 1.7;  % Standard deviation
factory_filt = imgaussfilt(factory, std, 'FilterSize', 9);

%% White Noise Addition

% Convert SNR to linear scale
snr_db = 10;
snr_lin = 10^(snr_db / 10);
image_pow = sum(factory_filt(:).^2 / numel(factory_filt));
noise_pow = image_pow / snr_lin;

factory_filt_noised = imnoise(factory_filt, 'gaussian', 0, noise_pow);

% figure, imshow(factory_filt_noised);
% [r, snr] = psnr(factory_filt_noised, factory);

%% Impulse Response and Frequency Response of Gaussian Kernel
impulse_response = imgaussfilt(impulse, std);

% Compute the frequency response using the FFT
frequency_response = fftshift(fft2(impulse_response));

% Compute the magnitude spectrum (absolute value)
magnitude_spectrum = abs(frequency_response);

show_img(factory, factory_filt_noised, impulse_response, magnitude_spectrum);


%% Wiener Denoise and Inverse Filter

gaussian_kernel = fspecial('gaussian', 9, std);

denoised_image = wiener2(factory_filt_noised, [5 5], noise_pow);
factory_wiener_inverse_restored = deconvwnr(denoised_image, gaussian_kernel);
% factory_wiener_inverse_restored = rescale(factory_wiener_inverse_restored);

% Apply a threshold to the filtered image
maxx = max(abs(factory_wiener_inverse_restored(:)));
threshold = 0.2 * maxx;
thresholded_image = factory_wiener_inverse_restored;
old_image = thresholded_image;
thresholded_image(abs(thresholded_image) < threshold) = 255;

figure, montage({factory_wiener_inverse_restored, thresholded_image}, 'BorderSize', 12)


%% Wiener Deconvolution

signal_var = var(factory(:));
nsr = noise_pow / signal_var;

factory_wiener_deconv = deconvwnr(factory_filt_noised, gaussian_kernel, nsr);

%% With Estimations

% Extract a representative noise patch from the image
noisePatch = factory_filt_noised(100:200, 100:200);

% Calculate the power of the noise patch
noisePower = var(noisePatch(:));

% Estimate the power of the signal
signalPower = mean(factory(:)).^2;

% Calculate the Noise Power Ratio (NPR)
NPR = noisePower / signalPower;

% Determine the Wiener filter parameters
factory_wiener_deconv_estimation = deconvwnr(factory_filt_noised, gaussian_kernel, NPR);

%%
figure, montage({factory_filt_noised, factory_wiener_deconv, factory_wiener_deconv_estimation})
%%
function show_img(image, image_filtered, impulse, freq_resp)
    figure;
    subplot(3, 2, 1);
    imshow(image);
    title('Original Image');
    
    subplot(3, 2, 2);
    imshow(image_filtered);
    title('Filtered Image');
    
    subplot(3, 2, 3);
    surf(impulse);
    title('Kernel Impulse Response');
    
    subplot(3, 2, 4);
    surf(freq_resp);
    title('Kernel Frequency Response');
    
%     Adjust the subplot arrangement and spacing
%   set(gcf, 'Position', [0.2, 0.2, 0.6, 0.6]);
end