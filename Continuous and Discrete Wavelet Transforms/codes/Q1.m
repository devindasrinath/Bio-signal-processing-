clearvars; close all; clc;

fs = 250;           % Sample frequency
N = 3000;           % Data length
t = (-N:N)/fs;      % Time scale
s = 0.01:0.1:2;     % Values of scaling factor (20 scales)


%% Maxican hat wavelet
num_scaling_factors = length(s);

% lets tau = 0
figure;
for i = 1:num_scaling_factors
    %generate wavelets
    wavelt(:, i) = (2/(sqrt(3*s(i))*(pi^(1/4))))*(1-(t/s(i)).^2).*exp(-((t/s(i)).^2) /2);
    
    %plot wavelets
    subplot(5,4,i);
    plot(t, wavelt(:, i));
    title(['Scale = ', num2str(s(i))]);
    xlabel('Time(s)'), ylabel('Amplitude');
    axis([-5 5 -1.5 3]);
    
    %calculate mean and energy values
    syms t_i;
    wavelt_M(i) = int((2/(sqrt(3*s(i))*(pi^(1/4))))*(1-(t_i/s(i)).^2).*exp(-((t_i/s(i)).^2) /2), 't_i', -inf, inf); %for mean
    wavelt_E(i) =  int(((2/(sqrt(3*s(i))*(pi^(1/4))))*(1-(t_i/s(i)).^2).*exp(-((t_i/s(i)).^2) /2))^2, 't_i', -inf, inf); %for energy
end


% ploting mean and energy
figure;
scatter(s, wavelt_M,'r');
hold on;
scatter(s, wavelt_E,'g');
hold off;
axis([0 2 -0.5 1.5]);
title('Mean and Energy with scaling factor');
legend('Mean', 'Energy'), xlabel('Scale Factor'), ylabel('Amplitude');


% Generating Spectra of wavelets
figure;
for i = 1:num_scaling_factors
    %generate FFT
    Fwavelt = fft(wavelt(:, i))/length(wavelt(:, i));
    hz = linspace(0,fs/2,floor(length(wavelt(:, i))/2)+1);

    %plot FFT
    subplot(5,4,i);
    plot(hz,2*abs(Fwavelt(1:length(hz))))
    title(['Scale = ', num2str(s(i))]), xlabel('Frequency (Hz)'), ylabel('Amplitude');
    %for i ==1 frequecncy distributed in ver large range. but in other
    %cases range is ve low.
    if (i==1)
        axis([0 100 0 0.15]);
    elseif(i>1 && i<=4)
        axis([0 10 0 0.15]);
    elseif(i>4)
        axis([0 5 0 0.15]);
    end
end



