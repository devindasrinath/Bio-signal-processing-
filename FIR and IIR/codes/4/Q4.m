%==================Designing FIR filters using windows=====================
%==========================================================================
%==========================================================================

clear all; clc;

%=========== Characteristics of window functions (use the fdatool)=========

%design Rectangular window for (M = 5,50,100)
filterDesigner;

%create windows wchi size is 51 , order is 50
win_rect = rectwin(51);
win_hann = hann(51);
win_hamm = hamming(51);
win_black = blackman(51);

%plot window functions
figure;
plot(0:50 , win_hann,0:50 , win_rect, 0:50 , win_hamm,0:50 , win_black);
title('Window functions');
ylim([0,1.2]);
legend('Rectangular' , 'Hanning' ,'Hamming' , 'Blackman');
xlabel('samples(n)');
ylabel('Amplitude');

%answers for otheres quastion in 4.1 part is available in .fda files in
%this folder.

%=========FIR Filter design and application using the Kaiser window========

fs = 500;
load ECG_with_noise.mat;

%Plot ECG signal
x = 1 : length(nECG);  %X axis
t = x / fs; %scaled X axis according to fs to obtain real time axis
figure;
plot(t,nECG);  %plot ECG signal
xlabel('Time(s)');
ylabel('ECG amplitude(mV)');
title('ECG signal with noise in time domain');


%=================Plot the power spectral density (PSD) estimate===========
figure
periodogram(nECG,[],length(nECG),fs);
title('Periodogram PSD estimate of nECG');



%===================Kaiser window for low pass =============================
f_pass = 123; %2.5Hz
f_stop = 127; %7.5Hz
delta = 0.001;

%calculate low pass filter paramerts
wp = 2*pi*(f_pass/fs);
ws = 2*pi*(f_stop/fs);

A = -20*log10(delta);
dw = ws - wp;
N_temp = ceil(((A-8)/(2.285*dw)) + 1);
N = (rem(N_temp,2) == 0) + N_temp;  %if length is even we should add 1
wc=(ws+wp)/2;
M = N -1;

disp( ['low pass order :', num2str(M)]);

% since A = 80 (A>50) , beta = 0.1102(A-8.7) is used
beta = 0.1102*(A-8.7);

Ib = besseli(0, beta);  % zeroth order modified Bessel function of the first kind

for n = 1:M+1
    % Calculating coefficients Kaiser window w(n)
    x = beta*sqrt(1-(((n-1)-M/2)/(M/2))^2);
    I0 = besseli(0,x);
    w(n) = I0/Ib;
    
    % Calculating coefficients of desired impulse response hd(n)
    if ((n-1)- (M/2)) == 0
        hd(n) = wc/pi;
    else 
        hd(n) = sin(wc*((n-1)-M/2))/(pi*((n-1)-M/2));
    end
end


%calculate coff of actual impulse response hd(n)
h_low = hd.*w;

gd_low = round(grpdelay(h_low,1 ,1));

%========FIR Filter design and application using the Kaiser window=========
f_pass = 7; %127.5Hz
f_stop = 3; %122.5Hz
delta = 0.001;

%calculate high pass filter paramerts
wp = 2*pi*(f_pass/fs);
ws = 2*pi*(f_stop/fs);

A = -20*log10(delta);
dw = - ws + wp;
N_temp = ceil(((A-8)/(2.285*dw)) + 1);
N = (rem(N_temp,2) == 0) + N_temp;  %if length is even we should add 1
wc=(ws+wp)/2;
M = N -1;

disp( ['High pass oreder :', num2str(M)]);

% since A = 80 (A>50) , beta = 0.1102(A-8.7) is used
beta = 0.1102*(A-8.7);

Ib = besseli(0, beta);  % zeroth order modified Bessel function of the first kind

for n = 1:M+1
    % Calculating coefficients Kaiser window w(n)
    x = beta*sqrt(1-(((n-1)-M/2)/(M/2))^2);
    I0 = besseli(0,x);
    w(n) = I0/Ib;
    
    % Calculating coefficients of desired impulse response hd(n)
    if ((n-1)- (M/2)) == 0
        hd(n) = (pi - wc)/pi;
    else 
        hd(n) = -1*sin(wc*((n-1)-M/2))/(pi*((n-1)-M/2));
    end
   
end


%calculate coff of actual impulse response hd(n)
h_high = hd.*w;

%visulaize the filter
fvtool(h_high,1 ,h_low, 1);
gd_high= round(grpdelay(h_high,1 ,1));

%==========apply filters to signals==================

y_temp = filter(h_low,1,nECG);
y_temp(1 : length(y_temp)- round(gd_low)) = y_temp(round(gd_low) + 1  : end);
y_temp(length(y_temp)- round(gd_low) + 1: end)=0;

%plot low pass filtered signal
figure;
plot(t,nECG,t,y_temp);  %plot ECG signal
xlabel('Time(s)');
ylabel('ECG amplitude(mV)');
title('Low pass Filtered ECG signal in time domain');
legend('nECG' ,'LPF');

%plot low pass filtered signal PSD
figure;
periodogram(y_temp,[],length(y_temp),fs);
title('Periodogram PSD estimate of low filtered ECG');


y_filtered = filter(h_high,1,y_temp);
y_filtered(1 : length(y_filtered)- round(gd_high)) = y_filtered(round(gd_high) + 1  : end);
y_filtered(length(y_filtered)- round(gd_high) + 1: end)=0;

%plot low ,high filtered signal
figure;
plot(t,nECG,t,y_filtered);  %plot ECG signal
xlabel('Time(s)');
ylabel('ECG amplitude(mV)');
title(' Low and High Filtered ECG signal in time domain');
legend('nECG' ,'LPF + HPF');

figure
periodogram(y_filtered,[],length(y_filtered),fs);
title('Periodogram PSD estimate of low,high filtered ECG');

%comobo filter design
f_comb1 = 50; %50Hz
f_comb2 = 100; %50Hz
f_comb3 = 150; %50Hz

w_comb1 = (f_comb1*2*pi)/fs;
w_comb2 = (f_comb2*2*pi)/fs;
w_comb3 = (f_comb3*2*pi)/fs;

%caluclate zeros
z_1 = cos(w_comb1) + 1i*sin(w_comb1);
z_1_c = cos(w_comb1) - 1i*sin(w_comb1);
z_2 = cos(w_comb2) + 1i*sin(w_comb2);
z_2_c = cos(w_comb2) - 1i*sin(w_comb2);
z_3 = cos(w_comb3) + 1i*sin(w_comb3);
z_3_c = cos(w_comb3) - 1i*sin(w_comb3);

%calculate polynomial coefficents
z_1_conv = conv( [1 -z_1] , [1 -z_1_c]);
z_2_conv = conv( [1 -z_2] , [1 -z_2_c]);
z_3_conv = conv( [1 -z_3] , [1 -z_3_c]);
z_12_conv = conv( z_1_conv , z_2_conv);
z_123_conv = conv( z_12_conv , z_3_conv);
z_123_conv_flip=flip(z_123_conv , 7);
 
%calculate polynomial gain
G = 1/polyval(z_123_conv,1);

comb = z_123_conv*G;

%visulaize comb filter
fvtool(comb,1);
y_filtered_notch = filter(comb , 1 ,y_filtered);

%plot filtered signal
figure;
plot(t ,nECG , t,y_filtered_notch);  %plot ECG signal
xlabel('Time(s)');
ylabel('ECG amplitude(mV)');
title(' Low, High and comb Filtered ECG signal in time domain');
legend('nECG' ,'LPF + HPF + comb');

%=================Plot the power spectral density (PSD) estimate===========
figure;

%extract PSD data of noisy and filter output
[pxx1,w1] = periodogram(nECG,[],length(nECG),fs);
[pxx2,w2] = periodogram(y_filtered_notch,[],length(y_filtered_notch),fs);

%plot PSD of noisy and filtered signal
plot(w1,10*log10(pxx1) , 'r',w2,10*log10(pxx2) , 'black')
title('Periodogram PSD estimate of nECG and filtered ');
xlabel('frequency(Hz)');
ylabel('Magnitude(dB)');
legend('nECG', 'Filtered ECG');

%visulaize combinede filter
cascade  = conv(conv(h_low ,h_high), comb);
fvtool(cascade,1);

