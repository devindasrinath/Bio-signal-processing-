%======================= FIR derivative filters============================
%==========================================================================
%==========================================================================

%=========FIR derivative filter properties (use the fvtool(b,a))===========
%==========================================================================
clear all; clc;

n_coff1 = [1,-1];
d_coff1 = 1;
n_coff2 = [1,0,-1];
d_coff2 = 1;
fvtool(n_coff1 ,d_coff1,n_coff2 ,d_coff2);


G1 = 0.5;
G2 = 0.5;
fvtool(n_coff1*G1 ,d_coff1,n_coff2*G2 ,d_coff2);

%==================FIR derivative filter application=======================
%==========================================================================

load  ECG_rec.mat; 

fs=128;
t = (0 : length(ECG_rec) - 1)/fs;

%adding noise to the raw ECG
snr =10; %10dB
nECG = awgn(ECG_rec,snr,'measured') + 2*sin((2*pi/4)*t) + 3*sin(((2*pi/2)*t) + (pi/4)) ;

%plot noisy ECG and raw ECG
figure;
plot(t ,nECG , 'black' , t ,ECG_rec,'r' );
xlim([0,15]);
xlabel('Time (s)');
ylabel('Voltage(mV)');
title('Noisy ECG waveform');
legend('nECG','raw ECg');

%first order derviatiove filter
y1 = filter(n_coff1,d_coff1,nECG);

%3 point central diff. derviative filter
y2 = filter(n_coff2,d_coff2,nECG);

%plot filtered ECG signals
figure;
plot(t ,ECG_rec ,'black',t ,y1 ,'b',t,y2,'r');
xlim([10,12]);
xlabel('Time (s)');
ylabel('Voltage(mV)');
title('raw ECG and filter outputs');
legend('ECG', 'First order','3 point central diff.');
