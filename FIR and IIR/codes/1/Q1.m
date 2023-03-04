%=======================Moving average MA(N) filter========================
%==========================================================================
%==========================================================================

%===========================Preliminaries==================================
%==========================================================================

%=========================Plot ECG signal==================================
fs = 500;
Data = load('ECG_template.mat');
ECG = Data.ECG_template;
x = 1 : length(ECG);  %X axis
t = x / fs; %scaled X axis according to fs to obtain real time axis
figure;
plot(t,ECG,'b');  %plot ECG signal
title('Typical ECG signal');
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');

%============Generate and plot nECG signal on same plot====================
snr = 5; %5dB
nECG = awgn(ECG,snr,'measured');
figure;
plot(t,ECG,'black',t,nECG,'r');  %plot ECG signal
title('5dB Gaussian noise added ECG signal');
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');
legend('ECG signal','ECG signal with noise');

%=================Plot the power spectral density (PSD) estimate===========
[pxx_ECG,f_ECG] = periodogram(ECG,[],2^nextpow2(length(ECG)),fs);
[pxx_nECG,f_nECG] = periodogram(nECG,[],2^nextpow2(length(nECG)),fs);
figure;
plot(f_ECG,10*log10(pxx_ECG) , 'b',f_nECG,10*log10(pxx_nECG) , 'r' )
title('PSD estimate of noisy ECG and ideal ECG');
xlabel('Frequecy(Hz)'), ylabel('Power/frequency(dB/Hz)');
legend('ECG signal','ECG signal with noise');


%========MA(3) filter implementation with a customised script==============
%==========================================================================

%=========================implement filter MA(3)===========================
N_3 = 3; %filter order
y_in = nECG;
y_in_padded = padarray(y_in,[0 N_3-1],0,'pre'); %pad left sides 0 
ma3ECG_1 = zeros(1,length(y_in)); 
for i = N_3 : length(y_in_padded)  %loop through all the indexes except padded ones
    sum = 0;
    for j = i : -1: i-(N_3-1)
        sum = sum + y_in_padded(j); %sum of N number of values
    end  
    ma3ECG_1(i - (N_3-1)) =  sum/3;   %averaging and save
end

%group delay = (order-1)/2 (for derived FIR filter)
%Therefore group delay = (3-1)/2 = 1
g_delay = (N_3-1)/2;

%=======compnasate group delay of filtered nECG signal=====================
% we can shift data or shift time axis for delay compansation proces.
% ma3ECG_1_t = zeros(1,length(ma3ECG_1));
% ma3ECG_1_t(1 : length(ma3ECG_1)- round(g_delay)) = ma3ECG_1(round(g_delay) + 1  : end);
% ma3ECG_1 = ma3ECG_1_t;
t_comp_ma3ECG_1 = t - g_delay/fs;

%===============Plot nECG compansated , filtered nECG signal===============
figure;
plot(t,nECG ,'b', t_comp_ma3ECG_1,ma3ECG_1 , 'r');  %plot nECG and ma3ECG_1 signal
xlabel('Time(s)'), ylabel('nECG amplitude(mV)');
title('nECG vs ma3ECG_1(compansated group delay)');
legend('nECG' , 'ma3ECG_1(compansated group delay)');
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');
grid on

%====Plot the power spectral density (PSD) estimate of nECG and ma3ECG_1===
figure;
xlabel('frequency(Hz)'), ylabel('Power(dB/Hz)');
[pxx_ma3ECG_1,f_ma3ECG_1] = periodogram(ma3ECG_1,[],length(nECG),fs);
plot(f_nECG,10*log10(pxx_nECG) ,'b', f_ma3ECG_1,10*log10(pxx_ma3ECG_1) , 'g');  %plot nECG and ma3ECG_1 signal
title('PSD estimate of nECG and ma3ECG_1');
legend('nECG' , 'ma3ECG_1(compansated group delay)');
xlabel('Frequecy(Hz)'), ylabel('Power/frequency(dB/Hz)');
grid on

%========MA(3) filter implementation using builtin function================
%==========================================================================

%==========================inbuilt filter MA(3)============================
N_3 = 3; %filter order
y_in = nECG;
n_coff = (1/N_3)*ones(1,N_3);
d_coff=1;
ma3ECG_2 = filter(n_coff,1,y_in);

%group delay = (order-1)/2 (for FIR filter)
%Therefore group delay = (3-1)/2 = 1
g_delay_3 = (N_3-1)/2;

%============compansate group delay of filtered nECG signal================
% ma3ECG_2_t(1 : length(ma3ECG_2)- round(g_delay)) = ma3ECG_2(round(g_delay) + 1  : end);
% ma3ECG_2_t(length(ma3ECG_2)- round(g_delay) + 1: end)=0;
% ma3ECG_1 = ma3ECG_1_t;
t_comp_ma3ECG_2 = t - g_delay_3/fs;

%=======Plot ECG ,nECG compansated and filtered nECG signal================
figure;
plot(t,ECG , 'b', t,nECG ,'r', t_comp_ma3ECG_2 ,ma3ECG_2 , 'g');  %plot nECG and ma3ECG_1 signal
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');
title('ECG ,nECG and delay compansated ma3ECG_2');
legend('ECG' , 'nECG' , 'Delay caompansated ma3ECG_2');

%=============================Use the fvtool(b,a)==========================
fvtool(n_coff,d_coff);
% 
% 
% 
%========MA(10) filter implementation using builtin function==============
%========================================================================

N_10 = 10;
n_coff = (1/N_10)*ones(1,N_10);
d_coff =1;

%=============================Use the fvtool(b,a)========================
fvtool(n_coff,d_coff);

%==========================filter nECG using MA(10)========================
ma10ECG = filter(n_coff,d_coff,y_in);

%=========================calculate group delay============================
g_delay_10 = (N_10-1)/2 ;

%=======compnasate group delay of filtered nECG signal=====================
% ma10ECG(1 : length(ma10ECG)- round(g_delay)) = ma10ECG(round(g_delay) + 1  : end);
% ma10ECG(length(ma10ECG)- round(g_delay) + 1: end)=0;
t_comp_ma10ECG = t - round(g_delay_10)/fs;

%=======Plot ECG ,nECG compansated and filtered nECG signal================
figure;
plot(t,ECG , 'b', t,nECG ,'r', t_comp_ma3ECG_2, ma3ECG_2 ,'g', t_comp_ma10ECG ,ma10ECG , 'black');  %plot nECG and ma3ECG_1 signal
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');
title('ECG ,nECG,  ma3ECG_2 and ma10ECG ');
legend('ECG' , 'nECG' , 'Delay compansated ma3ECG_2','Delay compansated ma10ECG');
xlim([0 , 0.7]);

%====================caculate MSE and find optimum filter==================
%==========================================================================

iterations =80;
mse = zeros(1,iterations);

for N = 1 : iterations
    %==========================inbuilt filter MA(N)========================
    y_in = nECG;
    n_coff = (1/N)*ones(1,N);
    d_coff=1;
    maNECG = filter(n_coff,d_coff,y_in);
    
    g_delay = (N-1)/2;
    
    %each time we have to time compansate signla accoridng to the group delay%
    maNECG(1 : length(maNECG)- round(g_delay)) = maNECG(round(g_delay) + 1  : end);
    maNECG(length(maNECG)- round(g_delay) + 1: end)=0;
    
    %=========================MSE calculation==============================
    diff  = (ECG - maNECG );
    mse(N) = (diff*diff') / length(ECG)  ;
end

%===================plot MSE vs order======================================
figure;
plot((1:iterations) , mse);
title('Mean squrred error(MSE) vs Filter order(n)');
xlabel('Filter order(n)'), ylabel('magnitude of MSE');

[min_mse,opt_order] = min(mse);

disp(['optimum order = ' , num2str(opt_order)] );

tic
maOPECG = filter((1/opt_order)*ones(1,opt_order),1,nECG);
toc
g_delay = (opt_order-1)/2;
t_comp_opt_MA = t - round(g_delay)/fs;


%========================Savitzky-Golay SG(N,L) filter===================
%========================================================================
%========================================================================

%==========================Application of SG(N,L)========================
%========================================================================

%==========================inbuilt filter SG(N,L)========================
N_3 = 3; %filter order
f_length = 11;
sg310ECG = sgolayfilt(nECG,N_3,1+(f_length*2));

%group delay = order/2 (for FIR filter)
%Therefore group delay = (3-1)/2 = 1
%g_delay = (N_3-1)/2;

%=======compnasate group delay of filtered nECG signal=====================
%sg310ECG(1 : length(sg310ECG)- round(g_delay)) = sg310ECG(round(g_delay) + 1  : end);
%sg310ECG(length(sg310ECG)- round(g_delay) + 1: end)=0;
%t_comp = t - (round(g_delay)) /fs;

%=======Plot ECG ,nECG and filtered nECG signal===============%
figure;
plot(t,ECG , 'b', t,nECG ,'r', t,sg310ECG , 'g');  %plot nECG and sg310ECG signal
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');
title('ECG ,nECG and sg310ECG');
legend('ECG' , 'nECG' , 'sg310ECG');

%==========================Optimum SG(N,L) filter parameters===============
iterations_N =70;
iterations_L =70;
mse = NaN(iterations_N, iterations_L) ;
min_mse = 100;
for N = 1 : iterations_N
    for L = 1 : iterations_L
        if N< 2*L        
            sgNLECG = sgolayfilt(nECG,N, 2*L + 1);            
            %g_delay = (N-1)/2;
            %each time we have to time compansate signla accoridng to the group delay%
            %sgNLECG(1 : length(sgNLECG)- round(g_delay)) = sgNLECG(round(g_delay) + 1  : end);
            %sgNLECG(length(sgNLECG)- round(g_delay) + 1: end)=0;

            %since time compansated signal have several zeros at the end of the
            %filtered signal we have set zeros at the end of the orginal signal
            %ECG_t = ECG;
            %ECG_t(length(ECG_t)- round(g_delay) + 1: end)=0;
            
            %=========================MSE calculation=================================%
            diff  = (ECG - sgNLECG );
            mse(N,L) = immse(ECG, sgNLECG); %(diff*diff') / (length(ECG) - round(g_delay));
            if min_mse > mse(N,L)
                min_mse = mse(N,L);
                op_N = N;
                op_L = L;
            end
        end
    end
end

figure
surf(1: iterations_N ,1:iterations_L ,mse);
ylabel('filter order(N)'),xlabel('filter length(L)'),zlabel('MSE');

disp(['Optimum N = ' , num2str(op_N)]);
disp(['Optimum L = ', num2str(op_L)]);

%==========================generate optimum filter=========================
N = op_N; %filter order
w_length = op_L*2 +1;
tic
sgOPECG = sgolayfilt(nECG,N,w_length);
toc

%Therefore group delay = (op_N-1)/2 = 1
%g_delay = (N-1)/2;

%=======compnasate group delay of filtered nECG signal=====================
%sgOPECG(1 : length(sgOPECG)- round(g_delay)) = sgOPECG(round(g_delay) + 1  : end);
%sgOPECG(length(sgOPECG)- round(g_delay) + 1: end)=0;
%t_comp_opt_SG = t - round(g_delay)/fs;

%======================Plot ECG ,sg310ECG and sgOPECG signal===============
figure;
plot(t,ECG , 'b',t,nECG , 'y', t ,sg310ECG , 'g', t,sgOPECG , 'r');  %plot sgOPECG and sg310ECG signal
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');
title('ECG ,sg310ECG and sgOPECG');
legend('ECG' ,'nECG', 'sg310ECG','sgOPECG');

 
%=======Plot ECG ,sg310ECG and sgOPECG signal===============%
figure;
plot(t,ECG , 'b',t,nECG , 'y', t_comp_opt_MA ,maOPECG , 'g', t ,sgOPECG , 'r');  %plot sgOPECG and sg310ECG signal
xlabel('Time(s)'), ylabel('ECG amplitude(mV)');
title('ECG ,maOPECG and sgOPECG');
legend('ECG' ,'nECG', 'maOPECG','sgOPECG');









