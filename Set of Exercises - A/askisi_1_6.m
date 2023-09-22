% Import the image
f = imread("clock.jpg");
f = rgb2gray(f);
image = im2double(f);

% whos f

%% Sobel Kernels for Edge Detection
% for v = 0.01:0.01:0.1
    v = 0.07;
    image_sobel_vertical = edge(image, 'Sobel', v,'vertical');
    image_sobel_horizontal = edge(image, 'Sobel', v,'horizontal');
    image_sobel_both = edge(image, 'Sobel', v, 'both');

    figure, montage({image, image_sobel_vertical, image_sobel_horizontal, image_sobel_both})
    title("Original image, Vertical Sobel, Horizontal Sobel, Both Directions Sobel ");
% end

%% Hough transform line detection

% Apply the Hough transform to detect lines
[H, theta, rho] = hough(image_sobel_both, 'thetaRes', 0.1);

figure, imshow(imadjust(mat2gray(H)), 'XData', theta, 'YData', rho, 'Border', 'loose')
daspect auto
axis on
xlabel('\theta')
ylabel('\rho')

% Find the peaks in the Hough accumulator matrix
peaks = houghpeaks(H,8, 'Threshold', 0.3*max(H(:)));
hold on
plot(round(theta(peaks(:,2))), rho(peaks(:,1)), 'Marker', 's', 'MarkerSize', 14, 'color', 'g')
% Extract the lines based on the peaks
hold off
lines = houghlines(image_sobel_both, theta, rho, peaks, 'FillGap', 25);

% Display the original image with detected lines
figure, imshow(image); hold on

% Plot the detected lines
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'y');
end

% Add title and labels
title('Edge Detection with Hough Transform');
xlabel('X');
ylabel('Y');