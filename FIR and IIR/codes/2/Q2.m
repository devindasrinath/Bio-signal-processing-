%============================Signal with multiple measurements=============
%==========================================================================
%==========================================================================

%================================Preliminaries=============================
%==========================================================================
clear all; clc;
load ABR_rec.mat; 

%plot whole ECG recording with stimuli
figure;
plot(ABR_rec);
xlabel('samples(n)');
ylabel('Amplitude(mV)');
title('Train of stimuli and ABRs');
legend('Stimuli','ABR train');

%plot zoomed ECG recording with stimuli
figure;
plot(ABR_rec);
xlabel('samples(n)');
ylabel('Amplitude(mV)');
title('zoomed Train of stimuli and ABRs');
legend('Stimuli','ABR train');
xlim([1000000,1003000]);

thresh = find(ABR_rec(:,1)>50);

%extract data using threshold
j=1;
for i=1:length(thresh)-1
if thresh(i+1)-thresh(i)>1; stim_point(j,1)=thresh(i+1); j=j+1;
end
end

%create epochs
j = 0;
for i=1:length(stim_point) j = j + 1;
epochs(:,j) = ABR_rec((stim_point(i)-80:stim_point(i)+399),2);
end

ensmbl_avg = mean(epochs(:,(1:length(stim_point))),2);

%plot portion of filtered signal
figure;
plot((-80:399)/40,ensmbl_avg);
xlabel('Time (ms)');
ylabel('Voltage(uV)');
title(['Ensemble averaged ABR from ',num2str(length(epochs)),' epochs']);

%=============================Improvement of the SNR=======================
%==========================================================================

%caculate MSE
[N,M] = size(epochs);

MSE = zeros(1, M);
for K = 1 : M
    y_k = mean(epochs(:,1:K) , 2);
    diff = (ensmbl_avg - y_k);
    MSE(K) = sqrt((diff'*diff)/N);
    SNR(K) = snr(ensmbl_avg,ensmbl_avg-y_k);
end

%plot MSE against number of epochs
figure;
plot(1:M,MSE);
xlabel('Number of epochs(k)');
ylabel('MSE_k');
title('MSE against number of epochs');

%MSE(dB) against number of epochs
figure;
plot(1:M,20*log10(MSE));
xlabel('Number of epochs(k)');
ylabel('MSE_k (dB)' );
title('MSE(dB) against number of epochs');

%plot SNR variation vs number of epochs
k= 1 : 1: M;
ideal_SNR =  10*log10(k)+SNR(1);
figure;
plot(k,ideal_SNR,'r',k,SNR,'b');
xlabel('Number of epochs(k)');
ylabel('SNR(dB)' );
title('SNR variation vs number of epochs');

%====================== Signal with repetitive patterns==================
%========================================================================

%================Viewing the signal and addition Gaussian white noise====
load ECG_rec.mat; 

%plot raw ECG data
fs=128;
t = (0 : length(ECG_rec) - 1)/fs;
figure;
plot(t ,ECG_rec);
xlabel('Time (s)');
ylabel('Voltage(mV)');
title('ECG wave datatset');

%plot portion of ECG recording
figure
plot(t ,ECG_rec);
hold on
xlim([0.2,2]);
xlabel('Time (s)');
ylabel('Voltage(mV)');
title('Zoomed ECG wave dataset');
hold off

%template sample range
x_1= 134;
x_2= 233;
ECG_template = ECG_rec(x_1:x_2);

%plot extracted ECG template
figure
plot(1: length(ECG_template) ,ECG_template);
xlabel('Sample(n)');
ylabel('Voltage(mV)');
title('Extracted ECG waveform(PQRST)');

%generate noisy signal
snr_ = 5; %5dB
nECG = awgn(ECG_rec, snr_, 'measured');

%plot noisy signal
figure
plot(t ,nECG);
xlabel('Time (s)');
ylabel('Voltage(mV)');
title('Noisy ECG recording');

%========Segmenting ECG into separate epochs and ensemble averaging========
ECG_template_padded = padarray(ECG_template,[0 (length(nECG) -length(ECG_template))] ,0,'post');

%plot cross corilation values
[xcorr_vals,lags] = xcorr(nECG,ECG_template_padded,'normalized');
xcorr_vals_extracted = xcorr_vals((length(xcorr_vals) + 1)/2 : length(xcorr_vals));
figure
plot(lags/fs, xcorr_vals)
xlabel('Lags (s)');
ylabel('Normalized value');
title('Normalized cross-correlation values');

% highest corrilations and there samples
%1 -> 37  , 0.08750
%2 -> 133 , 0.09574 
%3 -> 233 , 0.08709 
%4 -> 329 , 0.09189
%5 -> 419 , 0.09559
%6 -> 508 , 0.09896
%7 -> 601 , 0.08430
%8 -> 711 , 0.09832
%9 -> 831 , 0.08340

%extract highest corilation points%
corr_th = 0.08;
epoch_indexes_temp = zeros(1,length(xcorr_vals_extracted));

i = 1;
k = 1;
while i <=length(xcorr_vals_extracted)
    if xcorr_vals_extracted(i) > corr_th
        max = xcorr_vals_extracted(i);
        epoch_indexes_temp(k) = i;
        for j = 1 : 2
            if xcorr_vals_extracted(i+j) > max
                max = xcorr_vals_extracted(i+j);
                epoch_indexes_temp(k) = i+j;
            end
        end
    k=k+1;
    i= i + 2;
    end
    i=i+1;
end

epoch_indexes = epoch_indexes_temp(1 : k-1);
epochs = zeros(length(ECG_template) ,k-1);

%extract epochs related to highest corilation points%
k_new = 0;
for i = 1:k-1
    if epoch_indexes(i) + length(ECG_template) - 1 < length(nECG)
        epochs(: , i ) = nECG(epoch_indexes(i) : epoch_indexes(i) + length(ECG_template) - 1);
        k_new= k_new+1;
    end
end

epochs = epochs(:, 1 : k_new );
num_epochs = k_new;

[N,M] = size(epochs);

%calculate SNR
%MSE = zeros(1, M);
SNR = zeros(1, M);

for K = 1 : M
    y_k = mean(epochs(:,1:K) , 2);
    SNR(K) = snr(ECG_template' , y_k -ECG_template');
    %diff = (ECG_template' - y_k);
    %MSE(K) = sqrt((diff'*diff)/N);
end

%plot SNR aginst k
figure;
plot(1:M,SNR);
xlabel('Number of pulses(k)');
ylabel('SNR(dB)');
title('SNR against k');
                
%caculate ensemble average of arbiitarary selected teo epochs.
ensembled_1 = mean(epochs(:,1:5) , 2);
ensembled_2 = mean(epochs(:,1:70) , 2);

%plot SNR aginst k
t = (0 : length(ECG_template) - 1)/fs;
figure;
plot(t ,ECG_template , 'black' , t ,epochs(:,1), 'r', t ,ensembled_1 , 'b' ,t ,ensembled_2 , 'g');
xlabel('Time (s)');
ylabel('Voltage(mV)');
title('ECG waveform');   
legend('ECG','nECG','ensemble(5)','ensemble(70)');
        
%plot points of pulse detection 
len = fs*4;
x_axis = linspace(1.1,len,len)/fs;

%plot portion of ECG recording
figure('Name', 'xcorr pulse detection')
subplot(3,1,1);
plot(x_axis, ECG_rec(1:len))
title('ECG raw recording'),xlabel('Time (s)'), ylabel('Amplitude (mV)')

%plot delay adjusted cross corilation values
subplot(3,1,2);
plot(x_axis, xcorr_vals(floor(length(xcorr_vals)/2)+1:floor(length(xcorr_vals)/2+len)))
title('Adjusted xcorr values'),xlabel('Lag (s)'), ylabel('Normal Cross Correlation (mV)')

%plot noisy ECG with maxmimum cross corr. points
subplot(3,1,3);
plot(x_axis, nECG(1:len)), hold on
plot(epoch_indexes_temp(epoch_indexes_temp < len)/fs,xcorr_vals(epoch_indexes_temp < len),'*')
hold off
title('Noisy ECG and Pulse starting points'),xlabel('Time (s)'), ylabel('Amplitude (mV)')





