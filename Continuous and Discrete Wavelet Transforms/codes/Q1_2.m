clearvars; close all; clc;

fs = 250;           % Sample frequency
N = 3000;           % Data length
t = (-N:N)/fs;      % Time scale
s = 0.01:0.01:2;     % Values of scaling factor

num_scaling_factors = length(s);

%% generating the waveform

%prepare time axises
n1 = 1:1:(3*N/2) - 1; %should be exclude 3N/2
n2 = (3*N/2):1:3*N ; %should be include 3N/2

%generate main signal using 2 componenets at different time intervals
x = [sin(0.5*pi*n1/fs) sin(1.5*pi*n2/fs)];
n = [n1 n2];

figure('Name','x[n]')
plot(n, x, 'r');
title('x[n]'),xlabel('n'),ylabel('Amplitude');

%% Continuous Wavelet Decomposition
% lets tau = 0
for i = 1:num_scaling_factors
    %generate wavelet 
    wavelt = (2/(sqrt(3*s(i))*(pi^(1/4))))*(1-(t/s(i)).^2).*exp(-((t/s(i)).^2) /2);
    
    %convolute with x[n]
    conv_x = conv(x, wavelt, 'same'); 
    
    %store coefficents
    cwt_coff(i,:) = conv_x;
end

% [~,max_coff_index] = max((cwt_coff'),[],2);
% max_coff_scalar = max_coff_index*0.01;
% Ploting the spectrogram
figure;
h = pcolor(n, s, cwt_coff);
set(h, 'EdgeColor', 'none');
colormap jet
xlabel('Time (s)')
ylabel('Scale')
title('Spectrogram for x[n]')
% hold on;
% plot(n , max_coff_scalar);
% hold off;


%plot FFT of scalllignfactor  = 1.01 wavelet
% figure;
% wavelt_s_1 = (2/(sqrt(3*0.34)*(pi^(1/4))))*(1-(t/0.34).^2).*exp(-((t/0.34).^2) /2);
% Fwavelt = fft(wavelt_s_1)/length(wavelt_s_1);
% hz = linspace(0,fs/2,floor(length(wavelt_s_1)/2)+1);
% plot(hz,2*abs(Fwavelt(1:length(hz))))
% title(['Scale = ', num2str(1.01)]), xlabel('Frequency (Hz)'), ylabel('Amplitude');
% xlim([0 1.5]);