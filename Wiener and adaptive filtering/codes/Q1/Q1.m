clear all; clc;

%% load the ideal ECG
load idealECG.mat;
y_i  = idealECG - mean(idealECG);
fs = 500;
y_i_len = length(y_i);
t = linspace(0, y_i_len-1, y_i_len)*(1/fs);

%% create noisy signal
SNR =10; %10dB for awg
F_1 = 50; %50Hz for sin(power line noise)
x = awgn(y_i ,SNR,'measured' ) +  0.2*sin(2*pi*F_1*t); %noise + ideal
n_awg_n_50 = x - y_i; %noise 

%% plot all ECG signals
figure;
plot(t ,y_i);
title("ECG without noise"), xlabel('Time (s)'), ylabel('Amplitude (mV)');
figure;
plot(t, x);
title("ECG with noise"), xlabel('Time (s)'), ylabel('Amplitude (mV)');
figure;
plot( t, y_i);
title("ECG without noise(zoomed)"), xlabel('Time (s)'), ylabel('Amplitude (mV)');
xlim([0, 1.5]);
figure;
plot(t, x);
title("ECG with noise(zoomed)"), xlabel('Time (s)'), ylabel('Amplitude (mV)');
xlim([0,1.5]);

%% ===================== PART 1===============================%

%132 - 223 - > 92 samples
single_beat_period = 132 : 223;
y_desired = y_i(single_beat_period);
t_desired = t(single_beat_period);
% 
%121 - 143 ->23 samples
x_n_iso_period = 121:143;
x_n_iso = x(x_n_iso_period);

%replicate 4 times
noise = [x_n_iso x_n_iso x_n_iso x_n_iso];

%% plotting the desired ECG signal and the noise signal

figure('Name','desire ECG Single Beat and Noise estimate')
subplot(1,3,1)
plot(t_desired,y_desired)
title('Desired ECG Beat'), xlabel('Time (samples)'), ylabel('Amplitude (mV)')

subplot(1,3,2)
plot(1 : length(x_n_iso) ,x_n_iso,'g')
title('iso-noise'), xlabel('Time (samples)'), ylabel('Amplitude (mV)');

subplot(1,3,3)
plot(1:length(noise),noise,'r')
title('Replicated iso-noise'), xlabel('Time (samples)'), ylabel('Amplitude (mV)');

%% arbitarary filter order , let M  =18
M= 12;
W = weiner_weight_vector(y_desired, noise, 12);
disp(W);
%fvtool(W,1);

%apply weiner filter
y_hat = weiner_filter(x,W);

%plot filtered output and noise
n = 1 : length(y_hat);
figure;
plot(t , y_i ,  t, x, t , y_hat );
title("Desired signal , filter out and noisy signal| order =12 ");
legend("Ideal ECG" , "Noisy ECG", "Weiner Filtered ECG(order = 12) " );
xlabel('Time (s)'), ylabel('Amplitude (mV)')
figure;
plot(t , y_i ,t, x , t , y_hat );
title("Desired signal , filter out and noisy signal| order =12 ");
legend("Idela ECG", "Noisy ECG"  , "Weiner Filtered ECG(order = 12) ");
xlabel('Time (s)'), ylabel('Amplitude (mV)')
xlim([1, 1.5]);

%% find optimum filter order

% Find Optimum filter order and the coefficients

order_range = 50;
mse_values = NaN(1 , order_range);
for M = 2 : order_range
    W = weiner_weight_vector(y_desired, noise, M);
    y_hat = weiner_filter(x(single_beat_period),W); 
    mse_values(M) = immse(y_hat, y_desired);
end

[min_mse,opt_M] = min(mse_values);
disp(opt_M);

%plot mse with filter orders
figure;
scatter(opt_M, min_mse)
hold on;
plot(1 : order_range , mse_values);
hold off;
title("MSE vs Weigner filter order");
xlabel('Filter Order'), ylabel('MSE');

%apply opt ordered weiner filter
opt_W = weiner_weight_vector(y_desired, noise, opt_M);
y_hat_opt = weiner_filter(x,opt_W);
disp(opt_W);
%plot opt orderd filtered output and noise
figure;
plot(t , y_i  ,t, x,t , y_hat_opt );
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG(optimum order)')
title('Weiner Filtering of noisy ECG using Optimum order')
xlabel('Time (s)'), ylabel('Voltage (mV)');
figure;
plot(t , y_i  ,t, x, t , y_hat_opt );
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG(optimum order)')
title('Weiner Filtering of noisy ECG using Optimum order')
xlabel('Time (s)'), ylabel('Voltage (mV)');
xlim([1, 1.5]);

%visualize opitmum ordered fitler
fvtool(opt_W,1);

[Px_x, F1] = periodogram(x, [], [], fs);
[Px_noise,F2_n] = periodogram(n_awg_n_50, [], [], fs);
[Px_y_i,F3_i] = periodogram(y_i, [], [], fs);
[Px_y_hat,F4_y_hat] = periodogram(y_hat_opt, [], [], fs);

figure('Name','PSD')
semilogy( F2_n, Px_noise, F3_i, Px_y_i,F1, Px_x, F4_y_hat, Px_y_hat , 'black');
legend('Noise','Desired Signal','Noisy signal','Optimum Wiener Filtered Signal')
title('Power Spectral Density'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');
ylim([1.0000e-10 , 1]);




%% ======================= part 2===================%
linear_sig = zeros(1,92);
for i = 1: length(single_beat_period)
    if i < 10
        linear_sig(i) = 0;
    elseif i < 13
        linear_sig(i) = (0.05106 - (-0.1089))* (i - 10)/(13 - 10);  
    elseif i < 17
        linear_sig(i) = (0.05106 - (-0.1089));
    elseif i < 20
        linear_sig(i) = -(0.05106 - (-0.1089))* (i - 20)/(20 -17 );   
    elseif i < 30
        linear_sig(i) = 0;
    elseif i < 32
        linear_sig(i) =((-0.1089) - 0)* (i-30)/(32 - 30);
    elseif i < 36
        linear_sig(i) =((1.581- (-0.3489))* (i-32)/(36 - 32)) + (-0.3489);
    elseif i < 39
        linear_sig(i) = (((-0.5289) - (1.581))* (i-36)/(39 - 36)) + 1.581; 
    elseif i < 41 
        linear_sig(i) = (0- (-0.5289))* (i-41)/(41 - 39); 
    elseif i <  57
        linear_sig(i) = 0;
    elseif i < 66
        linear_sig(i) = (0.2611 - 0)* (i-57)/(66 - 57);  
    elseif i < 69
        linear_sig(i) = 0.2611;
    elseif i < 78 
        linear_sig(i) = ((0 - 0.2611 )* (i - 69)/(78 - 69)) + 0.2611; 
    else
        linear_sig(i) = 0;
    end
end
linear_sig = linear_sig - ones(1,92)*0.012;
figure;
plot(1: length(single_beat_period), linear_sig);
title('Linear ECG single Beat model'), xlabel('Sample'), ylabel('Amplitude (mV)');


%132 - 223 - > 92 samples
single_beat_period = 132 : 223;
t_desired = t(single_beat_period);

%121 - 143 ->23 samples
x_n_iso_period = 121 : 143;
x_n_iso = x(x_n_iso_period);

%replicate 4 times
noise = [x_n_iso x_n_iso x_n_iso x_n_iso];

%% plotting the desired ECG signal and the noise signal

figure('Name','desire ECG Single Beat and Noise estimate')
subplot(1,3,1)
plot(t_desired,linear_sig)
title('Desired ECG Beat'), xlabel('Time (s)'), ylabel('Amplitude (mV)')

subplot(1,3,2)
plot(1 : length(x_n_iso) ,x_n_iso,'g')
title('iso-noise'), xlabel('Time (s)'), ylabel('Amplitude (mV)');

subplot(1,3,3)
plot(1:length(noise),noise,'r')
title('Replicated iso-noise'), xlabel('Time (s)'), ylabel('Amplitude (mV)');


%% arbitarary filter order , let M  =12
M = 12;
W = weiner_weight_vector(linear_sig, noise, 12);
disp(W);
%fvtool(W,1);

%apply weiner filter
y_hat = weiner_filter(x,W);

%plot filtered output and noise
n = 1 : length(y_hat);
figure;
plot(t , y_i ,  t, x, t , y_hat );
title("Desired signal , filter out and noisy signal(order =12) ");
legend("Ideal ECG" , "Noisy ECG", "Weiner Filtered ECG(order = 12,linear model based) " );
xlabel('Time (s)'), ylabel('Amplitude (mV)')
figure;
plot(t , y_i ,t, x , t , y_hat );
title("Desired signal , filter out and noisy signal| order =12 ");
legend("Idela ECG", "Noisy ECG"  , "Weiner Filtered ECG(order = 12,linear model based) ");
xlabel('Time (s)'), ylabel('Amplitude (mV)')
xlim([1, 1.5]);

%% find optimum filter order for part 2

% Find Optimum filter order and the coefficients

order_range = 92;
mse_values = NaN(1 , order_range);
for M = 2 : order_range
    W = weiner_weight_vector(linear_sig, noise, M);
    y_hat = weiner_filter(x(single_beat_period),W); 
    mse_values(M) = immse(y_hat, y_i(single_beat_period));
end

[min_mse,opt_M] = min(mse_values);
disp(opt_M);

%plot mse with filter orders
figure;
scatter(opt_M, min_mse)
hold on;
plot(1 : order_range , mse_values);
hold off;
title("MSE vs Weigner filter order");
xlabel('Filter Order'), ylabel('MSE');

%apply opt ordered weiner filter
opt_W = weiner_weight_vector(linear_sig, noise, opt_M);
y_hat_opt = weiner_filter(x,opt_W);
disp(opt_W);
%plot opt orderd filtered output and noise
figure;
plot(t , y_i  ,t, x,t , y_hat_opt );
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG(optimum order,linear model based)')
title('Weiner Filtering of noisy ECG using Optimum order')
xlabel('Time (s)'), ylabel('Voltage (mV)');
figure;
plot(t , y_i  ,t, x, t , y_hat_opt );
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG(optimum order,linear model based)')
title('Weiner Filtering of noisy ECG using Optimum order')
xlabel('Time (s)'), ylabel('Voltage (mV)');
xlim([1, 1.5]);

%visualize opitmum ordered fitler
fvtool(opt_W,1);

[Px_x, F1] = periodogram(x, [], [], fs);
[Px_noise,F2_n] = periodogram(n_awg_n_50, [], [], fs);
[Px_y_i,F3_i] = periodogram(y_i, [], [], fs);
[Px_y_hat,F4_y_hat] = periodogram(y_hat_opt, [], [], fs);


figure;
semilogy( F2_n, Px_noise, F3_i, Px_y_i,F1, Px_x, F4_y_hat, Px_y_hat , 'black');
legend('Noise','Desired Signal','Noisy signal','Optimum Wiener Filtered Signal')
title('Power Spectral Density(linear model)'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');
ylim([1.0000e-10 , 1]);



%% 1.2 Frequency domain implementation of the Wiener filter

[y_hat_f, W_f] = weiner_fil_freq(y_i, n_awg_n_50, x);
figure;
plot(t, y_i, t, x, t, y_hat_f)
legend('y_i ECG','Noisy ECG','Weiner Filtered ECG(Frequency domain)')
title('Weiner Filtering Freq. domain'),
xlabel('Time (s)'),ylabel('Voltage (mV)');
xlim([1, 1.5]);
%% Plotting the spectrum

[Px_x, F1_x] = periodogram(x, [], [], fs);
[Px_noise,F2_noise] = periodogram(n_awg_n_50, [], [], fs);
[Px_y_i,F3_y_i] = periodogram(y_i, [], [], fs);
[Px_y_hat_f, F4_y_hat_f] = periodogram(y_hat_f, [], [], fs);

figure('Name','PSD Weiner Filter Frequency domain')
semilogy(F1_x, Px_x, F2_noise, Px_noise, F3_y_i, Px_y_i, F4_y_hat_f, Px_y_hat_f);
legend('Noisy signal','Noise','Desired Signal','Frequency domain Wiener Filtered Signal')
title('Power Spectral Density'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');
ylim([1.0000e-10 , 1]);
%% Comparing Frequency domain and Time domain weiner filter

figure('Name', 'Comparing Frequency domain and Time domain weiner filter')
plot(t, y_i, 'g', t, y_hat_opt,'b', t, y_hat_f, 'r')
legend('Ideal ECG','Filtered by Optimum Time Domain derived Weiner ', 'Filtered by Freq Domain derived Weiner Filter')
title('Comparing Frequency domain and Time domain weiner filter')
xlabel('Time (s)'), ylabel('Voltage (mV)');
xlim([1, 1.5]);

mse_time = immse(y_hat_opt, y_i);
mse_freq = immse(y_hat_f, y_i);

disp('Mean Square Error (Time domain)');
disp(mse_time);
disp('Mean Square Error (Frequency domain)');
disp(mse_freq);


%% 1.3 Effect on non-stationary noise on the Wiener filtering

% creating non stationary noise

f1 = 50;
f2 = 100;

t_p1 = t(1 : floor(length(t)/2));
t_p2 = t(floor(length(t)/ 2) + 1 : end);

n50_p1 = 0.2*sin(2*pi*f1*t_p1);
n100_p1 = 0.3*sin(2*pi*f2*t_p2);

non_stat_noise = [n50_p1 n100_p1];
non_sta_x = awgn(y_i ,SNR,'measured' ) +  non_stat_noise; %noise + ideal

N = length(non_sta_x);  
S_Xf  = fft(non_sta_x, N*2-1);
S_Yhat = W_f.* S_Xf;                    % Signal estimate from observation and using Wiener filter
y_hat_time = ifft(S_Yhat);              % converting to time domain
y_hat_non_stat = y_hat_time(1 : N);


figure('Name','Effect of Non Stationary Noise Comparison')
plot(t, y_i, t, y_hat_non_stat, t, y_hat_f)
xlim([t(1),t(end)])
legend('Ideal ECG','Non-Stationary Noise - Filtered','Stationary Noise - Filtered')
title('Effect of Non Stationary Noise Comparison after filtering with Weiner Freq')
xlabel('Time (s)'), ylabel('voltage(mV)')
xlim([7, 8.5]);

noise = y_hat_non_stat - y_i;

[Px_x, F1_x] = periodogram(non_sta_x, [], [], fs);
[Px_noise,F2_noise] = periodogram(noise, [], [], fs);
[Px_y_i,F3_y_i] = periodogram(y_i, [], [], fs);
[Px_y_hat_f, F4_y_hat_f] = periodogram(y_hat_non_stat, [], [], fs);

figure('Name','PSD Weiner Filter Frequency domain ( With Non Sattionary Noise)')
semilogy(F1_x, Px_x, F2_noise, Px_noise, F3_y_i, Px_y_i, F4_y_hat_f, Px_y_hat_f);
title('Power Spectral Density (With Non Sattionary Noise)'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');
legend('Non stationary Noise corrupted signal','Non Sattionary Noise','Desired Signal','Frequency domain Wiener Filtered Signal')
ylim([1.0000e-10 , 1]);