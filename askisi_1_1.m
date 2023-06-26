% Import the image
f = imread("aerial.tiff");
whos f

% Display the image
figure, imshow(f,[]);

% imtool(f)
%% 1.1.1 Γραμμική και Λογαριθμική Απεικόνιση του πλάτους του 2D Fourier
% Fourier Transform
Fourier_transform = fft2(f);

% Fourier Spectrum
fourier_spectrum = abs(Fourier_transform);
figure, imshow(fourier_spectrum, [])

% Fourier Tranformed Centered
Fourier_transform_centered = fftshift(Fourier_transform);
fourier_spectrum_centered = abs(Fourier_transform_centered);

% Fourier Spectrum Centered
figure, imshow(fourier_spectrum_centered, [])

% Log transformation of Fourier Spectrum
fourier_spectrum_log = log(1 + fourier_spectrum_centered);
figure, imshow(fourier_spectrum_log, [])

% Phase Angle
phase_angle = angle(Fourier_transform);
figure, imshow(phase_angle, [])

%The Polar form can now be expressed with the phase angle and the spectrum

%% 1.1.2 Χρήση Κατωπερατού Φίλτρου με padding
% Convert Image to float
Image = im2double(f);
[M,N] = size(Image);

% Padding Parameters
Padding_Parameters = 2*size(Image);


% Filter Parameters
n = 1;
D0 = 50;

% Butterworth Filter Generation
% Find the Meshgrid Arrays
u = 0:(Padding_Parameters(1) - 1);
v = 0:(Padding_Parameters(2) - 1);
idx = find(u> Padding_Parameters(1)/2);
u(idx) = u(idx) - Padding_Parameters(1);
idy = find(v> Padding_Parameters(2)/2);
v(idy) = v(idy) - Padding_Parameters(2);
[V,U] = meshgrid(v,u);

% Find the euclidian distance of the two meshgrids
Distance = hypot(U,V);

Butterworth_lp_filter = 1./(1 + (Distance./D0).^(2*n));

% Padding the image to avoid wraparound error
image_padded = padarray(Image, [size(Butterworth_lp_filter,1) - M, size(Butterworth_lp_filter,2) - N], "replicate", 'post');
figure, imshow(image_padded);

% Fast Fourier Transform of padded image
Fourier_padded = fft2(image_padded);

% Multiplication of H() and F()
G_Function = Butterworth_lp_filter.*Fourier_padded;

% Inverse FFT to get spatial domain 
g_Image = ifft2(G_Function);

% Crop to the size of the original image
g_Image = g_Image(1:M,1:N);

figure, imshow(g_Image);


%% 1.1.2 Χρήση Ανωπερατού Φίλτρου με padding
% Padding Parameters
Padding_Parameters = 2*size(Image);


% Filter Parameters
n = 1;
D0 = 50;

% Butterworth Filter Generation

Butterworth_hp_filter = 1.0 - (1./(1 + (Distance./D0).^(2*n)));

% Multiplication of H() and F()
G_Function_hp = Butterworth_hp_filter.*Fourier_padded;

% Inverse FFT to get spatial domain 
g_image_hp = ifft2(G_Function_hp);

% Crop to the size of the original image
g_image_hp = g_image_hp(1:M,1:N);

% Display Image
figure, imshow(g_image_hp);

% Intensity Scaling
g_image_scaled = mat2gray(g_image_hp);
figure, imshow(g_image_scaled)


%% High Frequency Emphasis Filtering
if false
    High_emphasis_filter = 0.5 + 2.0*Butterworth_hp_filter;
    G_Function_hp_emphasis = High_emphasis_filter.*Fourier_padded;
    g_image_hp_emphasis = ifft2(G_Function_hp_emphasis);
    g_image_hp_emphasis = g_image_hp_emphasis(1:M,1:N);
    figure, imshow(g_image_hp_emphasis)
    g_image_emphasis_scaled = mat2gray(g_image_hp_emphasis);
    figure, imshow(g_image_emphasis_scaled)
end

%% Απόκριση Πλάτους και Κρουστική Απόκριση των φίλτρων
% Magnitude Response of low-pass Butterworth Filter
figure, mesh(abs(fftshift(Butterworth_lp_filter)))
title("Magnitude Response of low-pass Butterworth")

% Magnitude Response of high-pass Butterworth Filter
figure, mesh(abs(fftshift(Butterworth_hp_filter)))
title("Magnitude Response of high-pass Butterworth")

% Impulse response of filter
impulse_response_lp = ifft2(Butterworth_lp_filter);
impulse_response_hp = ifft2(Butterworth_hp_filter);

% Impulse Response of low-pass Butterworth Filter
figure, mesh(fftshift(impulse_response_lp))
title("Impulse Response of low-pass Butterworth")

% Impulse Response of low-pass Butterworth Filter
figure, mesh(fftshift(impulse_response_hp))
title("Impulse Response of high-pass Butterworth")

if false
    figure, freqz2(impulse_response_lp)
    figure, freqz2(impulse_response_hp)
end