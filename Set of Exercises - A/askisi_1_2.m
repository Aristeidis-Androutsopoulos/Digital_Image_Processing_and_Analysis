% Import the image
f = imread("lenna.jpg");
f = rgb2gray(f);

% Resize the image to 255x255 using imresize
image = imresize(f, [256, 256], 'bilinear');
image_double = im2double(image);
% image = im2double(image);
whos f
figure, imshow(image);
min_value = min(f(:));
max_value = max(f(:));

%% Initializations

mse_threshold = 0;
mse_zonal = 0;
mse_zonal_after = 0;

ssim_threshold = 0;
ssim_zonal = 0;
ssim_zonal_after = 0;

dct_sum = zeros(32, 32);
dct_var = zeros(32, 32);

[height, width] = size(image);

num_sections_x = floor(width / 32);
num_sections_y = floor(height / 32);

%% Zonal Global Mask based on Variance

for y = 1:num_sections_y
   for x = 1:num_sections_x
        % Calculate start and end indices
        start_y = (y - 1) * 32 + 1;
        end_y = y * 32;
        start_x = (x - 1) * 32 + 1;
        end_x = x * 32;
        % Perform DCT to image section
        img_section = image(start_y:end_y, start_x:end_x);
        img_section_dct = dct2(img_section);
        % Sum all dct values
        dct_sum = dct_sum + img_section_dct;
   end
end

dct_mean = dct_sum ./ (num_sections_x * num_sections_y);
dct_var = ((dct_sum - dct_mean).^2) ./ ((num_sections_x * num_sections_y)-1);

%% Compression Loop for 0.05 to 0.5

% Define the compression ratio
for v = 0.05:0.05:0.5
    ratio_p = v;
    
    % Mask based on biggest variance
    num_coeffs = numel(img_section_dct);
    var_sorted = sort(dct_var(:), 'descend');
    threshold_idx = ceil(ratio_p * num_coeffs);
    var_thresh = var_sorted(threshold_idx);

    mask = dct_var; mask(mask < var_thresh) = 0; mask(mask ~= 0) = 1;

    % Define the block size for the blockproc block evaluation
    block_size = [32, 32]; % Adjust as desired

    % Define the processing function to apply the DCT and compression
    block_fun_mask = @(block) compressDCT_zonal_mask_after(block, mask);
    block_fun_threshold = @(block) compressDCT_threshold(block, ratio_p);
    block_fun_zonal = @(block) compressDCT_zonal(block, ratio_p);

    % Apply block processing to the image using DCT and thresholding
    compressed_image_zonal_after = blockproc(image, block_size, block_fun_mask);
    compressed_image_threshold = blockproc(image, block_size, block_fun_threshold);
    compressed_image_zonal = blockproc(image, block_size, block_fun_zonal);

    % Display the compressed image
    idct_image_threshold = rescale(compressed_image_threshold);
    figure, imshow(idct_image_threshold);
    title(["Threshold", ratio_p]);

    % Display the compressed image
    idct_image_zonal = rescale(compressed_image_zonal);
    figure, imshow(idct_image_zonal);
    title(["Zonal", ratio_p]);
    
    % Display the compressed image
    idct_image_zonal_after = rescale(compressed_image_zonal_after);
    figure, imshow(idct_image_zonal_after);
    title(["Zonal After", ratio_p]);
    
    % MSE of the three images compared to the original
    Thres = im2uint8(idct_image_threshold);
    Zonal = im2uint8(idct_image_zonal);
    Zonal_after = im2uint8(idct_image_zonal_after);
    
    err_threshold = immse(Thres, image);
    err_zonal = immse(image, Zonal);
    err_zonal_after = immse(image, Zonal_after);
  
    mse_threshold = [mse_threshold err_threshold];
    mse_zonal = [mse_zonal err_zonal];
    mse_zonal_after = [mse_zonal_after err_zonal_after];
    
    
    ssimval_thres = ssim(Thres,image);
    ssimval_zonal = ssim(Zonal,image);
    ssimval_zonal_after = ssim(Zonal_after,image);
    
    ssim_threshold = [ssim_threshold ssimval_thres];
    ssim_zonal = [ssim_zonal ssimval_zonal];
    ssim_zonal_after = [ssim_zonal_after ssimval_zonal_after];
end

% Print MSE
figure;
hold on
x_axis = linspace(0.05, 0.5, 11);
plot(x_axis, mse_zonal(2:end), 'DisplayName','Zonal')
hold on
plot(x_axis, mse_threshold(2:end), 'DisplayName','Threshold')
hold on
plot(x_axis, mse_zonal_after(2:end), 'DisplayName','Zonal after')
title("MSE of DCT Compression")
xlabel("Percentage of kept frequencies")
ylabel("MSE")
legend('Location', 'northeast')

figure;
hold on
x_axis = linspace(0.05, 0.5, 10);
plot(x_axis, ssim_threshold(2:end), 'DisplayName','SSIM Threshold')
hold on
plot(x_axis, ssim_zonal(2:end), 'DisplayName','SSIM Zonal')
hold on
plot(x_axis, ssim_zonal_after(2:end), 'DisplayName','SSIM Zonal after')
title("SSIM of DCT Compression")
xlabel("Percentage of kept frequencies")
ylabel("SSIM")
legend('Location', 'southeast')

%% Compression block functions
function compressed_block = compressDCT_threshold(block, compression_ratio)
    
    % Apply DCT to the block
    dct_block = dct2(block.data);
    
    % Determine the number of coefficients to retain based on the compression ratio
    num_coefficients = numel(dct_block);
    num_retained = ceil(compression_ratio * num_coefficients);
    
    % Sort the DCT coefficients in descending order of magnitude
    [~, sorted_indices] = sort(abs(dct_block(:)), 'descend');
    
    % Set coefficients below the threshold to zero
    dct_block(sorted_indices(num_retained+1:end)) = 0;
    
    % Reconstruct the block using inverse DCT
    compressed_block = idct2(dct_block);
end

function compressed_block = compressDCT_zonal(block, compression_ratio)
    
    % Map the ratio to an integer of (32, -32) rounded 
    num = round(interp1([0.001,1],[32,-32],compression_ratio));
    
    % Apply DCT to the block
    dct_block = dct2(block.data);
    
    % Create the zonal mask
    A = ones(size(dct_block));
    zonal_mask = triu(A, num);
    zonal_mask = flip(zonal_mask,2);

    % Get the compressed return of the mask*block
    compressed_block = zonal_mask .* dct_block;  
    
    % Reconstruct the block using inverse DCT
    compressed_block = idct2(compressed_block);
    
end

function compressed_block = compressDCT_zonal_mask_after(block, mask)
    
    % Apply DCT to the block
    dct_block = dct2(block.data);
    
    % Get the compressed return of the mask*block
    compressed_block = mask .* dct_block;  
    
    % Reconstruct the block using inverse DCT
    compressed_block = idct2(compressed_block);
    
end


