clearvars; close all; clc;

fs = 512;
rng(1);
%% Create X1 signal

n1 = 0:1:512-1; %first 512 points
n2 = 512:1:1023; %last 512 points

x1 = [(2*sin(20*pi*n1/fs) + sin(80*pi*n1/fs))  (0.5*sin(40*pi*n2/fs) + sin(60*pi*n2/fs))];

figure;
plot([n1 n2], x1);
xlim([0 1023]);
title('signal x_1 [n]'), xlabel('Time(s)'), ylabel('Amplitude');

%%  Constructing X2 signal

x2 = zeros(1, 1024); %array with 1024 samples
for i = 0:1023
    if (i >= 0 &&  i<64)
        x2(i+1) = 1;
    elseif (i >= 192 &&  i<256)
        x2(i+1) = 2;
    elseif (i >= 256 &&  i<512)
        x2(i+1) = -1;
    elseif (i >= 512 &&  i<704)
        x2(i+1) = 3;
    elseif (i >= 704 &&  i<960)
        x2(i+1) = 1;
    else
        x2(i+1)=0;
    end
end

figure;
plot(x2);
axis([0 1024 -1.5 3.5]);
title('X_2 [n]'), xlabel('Time(s)'), ylabel('Amplitude');

%% Corrupt these signals with AWGN of 10 dB SNR

y1 = awgn(x1, 10,'measured');
y2 = awgn(x2, 10,'measured');

figure;
plot([n1 n2], x1, [n1 n2], y1);
legend('ideal signal x_1 [n]', 'Noisy Signal y_1 [n]'),
title('x_1 [n] , y_1 [n] (noisy)'), xlabel('Time(s)'), ylabel('Amplitude');
xlim([0 1024]);

figure;
plot([n1 n2], x2, [n1 n2], y2);
legend('ideal signal x_2 [n]', 'Noisy Signal y_2 [n]'),
title('x_2 [n] , y_2 [n] (noisy)'), xlabel('Time(s)'), ylabel('Amplitude');
axis([0 1024 -3.5 4.5]);



%% Observe the morphology of the wavelet and scaling functions of Haar and Daubechies tap 9 using wavefun() command and the waveletAnalyzer GUI.

% Harr Wavelet 
% basic usage  :
%      [phi, psi ,xval] = wavefun(name, iter);
%in this case lets iter  = 10
[phi_haar, psi_haar, xval_haar] = wavefun('haar', 10); 
figure;
subplot(1,2,1);
plot(xval_haar, psi_haar, 'black');
title('Wavelet Function (Haar)');
subplot(1,2,2);
plot(xval_haar, phi_haar, 'red');
title('Scaling Function (Haar)');

% Daubechies tap 9 Wavelet 
[phi_deb9,psi_deb9, xval_deb9] = wavefun('db9', 10); 
figure;
subplot(1,2,1);
plot(xval_deb9, psi_deb9,'green');
title('Wavelet Function (Daubechies tap 9)');
subplot(1,2,2);
plot(xval_deb9, phi_deb9, 'blue');
title('Scaling Function (Daubechies tap 9)');


%% 10 level Wavelet Decomposition

%% For  Signal Y1 

% case 1 : Haar Wavelet
%we use wavedec for decomposisiton
% basic usage  :
%      [c, l] = wavedec(signal, stpes, wavelt name);
%in this case level = 10
[c_haar_y1, l_haar_y1] = wavedec(y1, 10, 'haar');
A_harr_y1 = appcoef(c_haar_y1, l_haar_y1, 'haar'); %for approximation corifccents
[Dh1_1,Dh1_2,Dh1_3,Dh1_4,Dh1_5,Dh1_6,Dh1_7,Dh1_8,Dh1_9,Dh1_10] = detcoef(c_haar_y1, l_haar_y1, [1 2 3 4 5 6 7 8 9 10]);

%plot detail coeffcents
figure;
for i= 1:10
    D_harr_y1 = detcoef(c_haar_y1, l_haar_y1, i);
    subplot(11,1,i);
    stem(D_harr_y1,'Marker','.');
    title(['Level ' num2str(i) ' Decomposition of y_1 [n] using Haar wavelet']);
end

%plot approx. coeffcents
subplot(11,1,11);
stem(A_harr_y1);
title(['Level' num2str(10) ' Approximation Coefficients of y_1 [n]']);


% case 2 : Daubechies tap 9 Wavelet
% in this case level = 10
[c_db9_y1, l_db9_y1] = wavedec(y1, 10, 'db9');
A_db9_y1 = appcoef(c_db9_y1, l_db9_y1, 'db9');
[Ddb1_1,Ddb1_2,Ddb1_3,Ddb1_4,Ddb1_5,Ddb1_6,Ddb1_7,Ddb1_8,Ddb1_9,Ddb1_10] = detcoef(c_db9_y1, l_db9_y1, [1 2 3 4 5 6 7 8 9 10]);

%plot detail coeffcents
figure;
for i= 1:10
    D_db9_y1 = detcoef(c_db9_y1, l_db9_y1, i);
    subplot(11,1,i);
    stem(D_db9_y1,'Marker','.');
    title(['Level ' num2str(i) ' Decomposition of y_1 [n] using db9 wavelet']);
end

%plot approx. coeffcents
subplot(11,1,11);
stem(A_db9_y1);
title(['Level' num2str(10) 'Approximation Coefficients of y_1 [n]']);


%% For  Signal Y2
% case 1 : Haar Wavelet 

% case 1 : Haar Wavelet
%we use wavedec for decomposisiton
% basic usage  :
%      [c, l] = wavefun(signal, stpes, wavelt name);
%in this case level = 10
[c_haar_y2, l_haar_y2] = wavedec(y2, 10, 'haar'); 
A_harr_y2 = appcoef(c_haar_y2, l_haar_y2, 'haar'); %for approximation corifccents
[Dh2_1,Dh2_2,Dh2_3,Dh2_4,Dh2_5,Dh2_6,Dh2_7,Dh2_8,Dh2_9,Dh2_10] = detcoef(c_haar_y2, l_haar_y2, [1 2 3 4 5 6 7 8 9 10]);

%plot detail coeffcents
figure;
for i= 1:10
    D_harr_y2 = detcoef(c_haar_y2, l_haar_y2, i);
    subplot(11,1,i);
    stem(D_harr_y2,'Marker','.');
    title(['Level ' num2str(i) ' Decomposition of y_2 [n] using Haar wavelet']);
end

%plot approx. coeffcents
subplot(11,1,11);
stem(A_harr_y2);
title(['Level' num2str(10) ' Approximation Coefficients of y_2 [n]']);


% case 2 : Daubechies tap 9 Wavelet
% in this case level = 10
[c_db9_y2, l_db9_y2] = wavedec(y2, 10, 'db9');
A_db9_y2 = appcoef(c_db9_y2, l_db9_y2, 'db9');
[Ddb2_1,Ddb2_2,Ddb2_3,Ddb2_4,Ddb2_5,Ddb2_6,Ddb2_7,Ddb2_8,Ddb2_9,Ddb2_10] = detcoef(c_db9_y2, l_db9_y2, [1 2 3 4 5 6 7 8 9 10]);

%plot detail coeffcents
figure;
for i= 1:10
    D_db9_y2 = detcoef(c_db9_y2, l_db9_y2, i);
    subplot(11,1,i);
    stem(D_db9_y2,'Marker','.');
    title(['Level ' num2str(i) ' Decomposition of y_2 [n] using db9 wavelet']);
end

%plot approx. coeffcents
subplot(11,1,11);
stem(A_db9_y2);
title(['Level' num2str(10) ' Approximation Coefficients of y_2 [n]']);


%% Discrete Wave Resconstruction
% case y1 
% using haar waveltes
x1_haar_rec = wavereconst('haar',A_harr_y1,Dh1_1,Dh1_2,Dh1_3,Dh1_4,Dh1_5,Dh1_6,Dh1_7,Dh1_8,Dh1_9,Dh1_10);

figure;
plot(0:1023, x1_haar_rec);
title('y_1 [n] reconstructed using Haar wavelet');
xlabel('Samples(n)');
ylabel('Amplitude');

% using db9
x1_db9_rec = wavereconst('db9',A_db9_y1,Ddb1_1,Ddb1_2,Ddb1_3,Ddb1_4,Ddb1_5,Ddb1_6,Ddb1_7,Ddb1_8,Ddb1_9,Ddb1_10);

figure;
plot(0:1023, x1_db9_rec);
title('y_1 [n] reconstructed using db9 wavelet');
xlabel('Samples(n)');
ylabel('Amplitude');
% case X2 - haar
x2_haar_rec = wavereconst('haar',A_harr_y2,Dh2_1,Dh2_2,Dh2_3,Dh2_4,Dh2_5,Dh2_6,Dh2_7,Dh2_8,Dh2_9,Dh2_10);

figure;
plot(0:1023, x2_haar_rec);
title('y_2 [n] reconstructed using Haar wavelet');
xlabel('Samples(n)');
ylabel('Amplitude');
% X2 - db9
x2_db9_rec = wavereconst('db9',A_db9_y2,Ddb2_1,Ddb2_2,Ddb2_3,Ddb2_4,Ddb2_5,Ddb2_6,Ddb2_7,Ddb2_8,Ddb2_9,Ddb2_10);

figure;
plot(0:1023, x2_db9_rec);
title('y_2[n] reconstructed using db9 wavelet');
xlabel('Samples(n)');
ylabel('Amplitude');

%% Energy calculation to compare recontructed and original signal

E_y1 = sum(abs(y1).^2);
E_y1_haar = sum(abs(x1_haar_rec).^2);
E_y1_db9 = sum(abs(x1_db9_rec).^2);
E_y2 = sum(abs(y2).^2);
E_y2_haar = sum(abs(x2_haar_rec).^2);
E_y2_db9 = sum(abs(x2_db9_rec).^2);

disp(['E(y_1) = ', num2str(E_y1)]);
disp(['E(Reconstructed y_1 | haar wavelet) = ', num2str(E_y1_haar)]);
disp(['E(Reconstructed y_1 | db9 wavelet) = ', num2str(E_y1_db9)]);
disp (' ')
disp(['E(y_2) = ', num2str(E_y2)]);
disp(['E(Reconstructed y_2 | haar wavelet) = ', num2str(E_y2_haar)]);
disp(['E(Reconstructed y_2 | db9 wavelet) = ', num2str(E_y2_db9)]);

%%  Signal Denoising with DWT
% case : y1 - haar
% Thrhsold = 1
[x1_rec_harr ,c_sorted_decend_harr_1] = wavedenoise(  'haar' ,y1, 10, 1);  
plot_denoise_wave('harr' ,'x_1','y_1' , c_sorted_decend_harr_1, x1_rec_harr, x1);
rmse_x1_harr = sqrt(sum(abs((x1 - x1_rec_harr)).^2)/length((x1 - x1_rec_harr))); % Calculating the RMSE between the original and the reconstructed signal
disp(['RMSE x1 reconstructed | haar wavelet = ' num2str(rmse_x1_harr)]);


% case : y1 - db9
% Thrhsold = 1
[x1_rec_db9 ,c_sorted_decend_db9_1] = wavedenoise( 'db9' ,y1,  10, 1);  
plot_denoise_wave('db9' ,'x_1','y_1' , c_sorted_decend_db9_1, x1_rec_db9, x1);
rmse_x1_db9 = sqrt(sum(abs((x1 - x1_rec_db9)).^2)/length((x1 - x1_rec_db9))); % Calculating the RMSE between the original and the reconstructed signal
disp(['RMSE x1 reconstructed | db9 wavelet = ' num2str(rmse_x1_db9)]);

% case : y2 - haar
% Thrhsold = 2
[x2_rec_harr ,c_sorted_decend_harr_2] = wavedenoise( 'haar' , y2, 10, 2);  
plot_denoise_wave('harr' ,'x_2','y_1' , c_sorted_decend_harr_2, x2_rec_harr, x2);
rmse_x2_harr = sqrt(sum(abs((x2 - x2_rec_harr)).^2)/length((x2 - x2_rec_harr))); % Calculating the RMSE between the original and the reconstructed signal
disp(['RMSE x2 reconstructed | haar wavelet = ' num2str(rmse_x2_harr)]);


% case : y2 - db9
% Thrhsold = 2
[x2_rec_db9 ,c_sorted_decend_db9_2] = wavedenoise( 'db9' ,y2,10, 2);  
plot_denoise_wave('db9' ,'x_2','y_1', c_sorted_decend_db9_2, x2_rec_db9, x2);
rmse_x2_db9 = sqrt(sum(abs((x2 - x2_rec_db9)).^2)/length((x2 - x2_rec_db9))); % Calculating the RMSE between the original and the reconstructed signal
disp(['RMSE x2 reconstructed | db9 wavelet = ' num2str(rmse_x2_db9)]);

%% Signal Compression

%load data
load('ECGsig.mat');                    
fs = 257;                               % sampling freqency 
         
%caculate the maximum number of levels can be applied . in this case we use this level for proper decompossing
num_levels = ceil(log2(length(aVR)));
percentage = 0.99;

%plot loaded signal
figure;
plot(1:(length(aVR)), aVR);
title('aVR Lead of ECG signal (fs = 257 Hz)');
xlabel('Samples(n)'), ylabel('Voltage (mV)');
xlim([0 length(aVR)]);


% case 1: haar wavelet 

% Obtain wavelet coefficients of aVR
[c_haar, l_haar] = wavedec(aVR, num_levels, 'haar');
approx_haar = appcoef(c_haar, l_haar, 'haar');

%plot detailed and apporiximation oefficents
figure;
for i= 1:num_levels
    %extract detailed coefficents
    haar_D = detcoef(c_haar, l_haar, i);
    
    %plot detailed coefficents
    subplot(num_levels+1,1,i);
    stem(haar_D,'Marker','.');
    title(['Level ' num2str(i) ' Decomposition of aVR signal | Haar wavelet']);
end
subplot(num_levels+1,1,num_levels+1);
stem(approx_haar);
title(['Level' num2str(num_levels) ' Approximation Coefficients']);



%compress wavforms using haar wavelet
[x_compressed_haar ,c_sorted_decend_haar] = wavecompress('haar', aVR, num_levels, percentage);
plot_denoise_wave('haar' ,'aVR','aVR' , c_sorted_decend_haar, x_compressed_haar, aVR);
rmse_aVR_haar = sqrt(sum(abs((aVR - x_compressed_haar)).^2)/length(aVR)); % Calculating the RMSE between the original and the reconstructed signal
disp(['RMSE aVR reconstructed | haar wavelet = ' num2str(rmse_aVR_haar)]);

% case 2: db9 wavelet 

% Obtain wavelet coefficients of aVR
[c_db9, l_db9] = wavedec(aVR, num_levels, 'db9');
approx_db9 = appcoef(c_db9, l_db9, 'db9');

%plot detailed and apporiximation oefficents
figure;
for i= 1:num_levels
    %extract detailed coefficents
    db9_D = detcoef(c_db9, l_db9, i);
    
    %plot detailed coefficents
    subplot(num_levels+1,1,i);
    stem(db9_D,'Marker','.');
    title(['Level ' num2str(i) ' Decomposition of aVR signal | db9 wavelet']);
end
subplot(num_levels+1,1,num_levels+1);
stem(approx_db9);
title(['Level' num2str(num_levels) ' Approximation Coefficients'])

%compress wavforms using db9 wavelet
[x_compressed_db9 ,c_sorted_decend_db9] = wavecompress('db9', aVR, num_levels, percentage);
plot_denoise_wave('db9' ,'aVR' ,'aVR', c_sorted_decend_db9, x_compressed_db9, aVR);
rmse_aVR_db9 = sqrt(sum(abs((aVR - x_compressed_db9)).^2)/length(aVR)); % Calculating the RMSE between the original and the reconstructed signal
disp(['RMSE aVR reconstructed | db9 wavelet = ' num2str(rmse_aVR_db9)]);
