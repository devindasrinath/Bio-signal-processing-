clear all; clc;

%% generate signals

%generate main signal
%y_i , t , n , n_50 ,n_100 ,x are column vectors
N = 5000; % Number of points
t = linspace(0,5,N)'; % Time vector with fs = 500 Hz
y_i = sawtooth(2*pi*2*t(1:N,1),0.5); % Sawtooth signal

%generate noise componenets
n_50 = 0.2*sin(2*pi*50*t(1:N/2,1)); % Sinusoid at 50 Hz
n_100 = 0.3*sin(2* pi*100*t(N/2+1:N,1)); % Sinusoid at 100 Hz
snr = 10; %10dB
nwg = y_i - awgn(y_i,snr,'measured'); % Gaussian white noise
n = nwg + [n_50; n_100]; %n??(?) + ?50(?) + ?100(?)

%generate noisy signal
x = y_i + n; 

%generate referance signal
a = 1.987;
phi_1 = pi/3; %?1
phi_2 = pi/4; %?2
n_50_100_ph1_ph2 = [ (0.2*sin(2*pi*50*t(1:N/2,1) + phi_1)); (0.3*sin(2* pi*100*t(N/2+1:N,1) + phi_2)) ]; %?50(??50 + ?1) + ?100(??100 + ?2)
r = a*(nwg +  n_50_100_ph1_ph2);%?(???(?) + ?50(??50 + ?1) + ?100(??100 + ?2))

%plot all the input signals and expected outputs for  adaptive filtering
figure('Name', 'Signals use for Adaptive Filtering');
subplot(3,1,1);
plot(t, y_i);
title('Ideal signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(3,1,2);
plot(t, x);
title('Noisy signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(3,1,3);
plot(t, r);
title('Referance signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
linkaxes();



%% PART 1 (LMS Algorithm)
%find the best M and mu
M_max = 30;
M_min  = 1;
M_range = M_min : M_max;

mu_min = 0.001; %mu_max depend on the order;
mu_points = 100;
lambda_max= 20*M_max*((x'*x)/length(x));
mu_max = 2/lambda_max;
mu_range = linspace(mu_min, mu_max, mu_points);
mse_values = NaN(M_max , mu_points);
min_mse = 100;
for M = M_range 
    mu_index = 1;
    for mu = mu_range
        [y_hat, W] = LMS(x' , r', M, mu);
        mse_values(M,mu_index) = immse(y_hat', y_i);
        if(mse_values(M,mu_index) <min_mse)
            min_mse = mse_values(M,mu_index);
            opt_M1 = M;
            opt_mu = mu;
            y_hat_opt  = y_hat;
        end
        mu_index = mu_index+1;
        disp(['M = '  num2str(M) ' mu = '  num2str(mu)]);
        
    end
end

figure
surf( M_range,mu_range, mse_values');
title('MSE Variation when Adaptive Filtering(LMS) with diffrerent orders and mu values') 
colorbar
xlabel('M - Order'), ylabel('mu'),zlabel('MSE');
colormap('jet')

disp(['optimum M : ' num2str(opt_M1) ', optimum mu : ' num2str(opt_mu)]);


%plot all the input signals and expected outputs for  adaptive filtering
figure('Name', 'Signals use for Adaptive Filtering');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(4,1,1);
plot(t, y_i);
title('Ideal signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(4,1,2);
plot(t, x);
title('Noisy signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(4,1,3);
plot(t, y_hat_opt);
title('filtered signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(4,1,4)
plot(t,abs(y_hat_opt-y_i'));
title('Absolute Error');xlabel('Time (s)'), ylabel('Amplitude (mV)');
linkaxes();

%% PART 2 (RMS Algorithm)
%find the best M and lambda


M_max = 25;
M_min  = 1;
M_range = M_min : M_max;

lambda_min = 0.8; %mu_max depend on the order;
lambda_points = 100;
lambda_max=1;
lambda_range = linspace(lambda_min, lambda_max, lambda_points);

mse_values = NaN(M_max , lambda_points);
min_mse = 100;
opt_M = 0;
opt_lambda=0;
for M = M_range 
    lambda_index = 1;
    for lambda = lambda_range
        [y_hat, W] =  RLS(x' , r' , M,lambda);
        mse_values(M,lambda_index) = immse(y_hat', y_i);
        if(mse_values(M,lambda_index) <min_mse)
            min_mse = mse_values(M,lambda_index);
            opt_M2 = M;
            opt_lambda = lambda;
            y_hat_opt_2  = y_hat;
        end
        lambda_index = lambda_index+1;
        disp(['M = '  num2str(M) ' lambda = '  num2str(lambda)]);
        
    end
end

figure
surf( M_range,lambda_range, mse_values');
title('MSE Variation when Adaptive Filtering(LMS) with diffrerent orders and lambda values') 
colorbar
xlabel('M - Order'), ylabel('lambda'),zlabel('MSE');
colormap('jet')

disp(['optimum M : ' num2str(opt_M2) ', optimum lambda : ' num2str(opt_lambda)]);



%plot all the input signals and expected outputs for  adaptive filtering
figure('Name', 'Signals use for Adaptive Filtering');
subplot(4,1,1);
plot(t, y_i);
title('Ideal signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(4,1,2);
plot(t, x);
title('Noisy signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(4,1,3);
plot(t, y_hat_opt_2);
title('filtered signal');xlabel('Time (s)'), ylabel('Amplitude (mV)');
subplot(4,1,4)
plot(t,abs(y_hat_opt_2-y_i'));
title('Absolute Error');xlabel('Time (s)'), ylabel('Amplitude (mV)');
linkaxes();



%% Comparison between LMS_method and RLS_method at their best outputs

mu_cmp = 0.001;
M_cmp_lms = 15;
[err1, W] = LMS(x', r', M_cmp_lms, mu_cmp);

lamda_cmp = 0.999;          
M_cmp_rls = 15;
[err2, W] = RLS(x', r', M_cmp_rls, lamda_cmp);

figure
subplot(2,1,1)
plot(t, abs(err1 - y_i'));
title(['Absolute error when using LMS algorithm mu = ' num2str(mu_cmp) ' M = ' num2str(15)  ]);
xlabel('Time(s)');
ylabel('Voltage (mV)');
grid on
subplot(2,1,2)
plot(t, abs(err2 - y_i'))
title(['Filtered Signal of the ANC filter using the RLS algorithm lambda = ' num2str(lamda_cmp) ' M = ' num2str(15) ]),
xlabel('Time(s)');
ylabel('Voltage (mV)');
grid on



%% Adaptive FIltering use for ECG signal

load('idealECG.mat')
fs = 500;
y_i = idealECG - mean(idealECG);
N = length(y_i);                                  
t = linspace(0,N/fs,N);                  
                  
%generate noise components
n_50 = 0.2*sin(2*pi*50*t(1 : N/2));  %50Hz sin component            
n_100 = 0.3*sin(2* pi*100*t(N/2+1 : N));   %100Hz sin component              
nwg = y_i - awgn(y_i, 10,'measured');  % Gaussian white noise of 10dB
x = y_i + nwg + [n_50 n_100];         % noisy ECG


% Generating refenerace signal
a = 1.987;
phi_1 = pi/4; %?1
phi_2 = pi/3; %?2
n_50_100_ph1_ph2 = [ (0.2*sin(2*pi*50*t(1:N/2) + phi_1)) (0.3*sin(2* pi*100*t(N/2+1:N) + phi_2)) ]; %?50(??50 + ?1) + ?100(??100 + ?2)
r = a*(nwg +  n_50_100_ph1_ph2);%?(???(?) + ?50(??50 + ?1) + ?100(??100 + ?2))

% apply LMS alogrithm to ECG_sig
[y_hat_lms_ecg, W] = LMS(x, r, 15, 0.001);

% apply RLS alogrithm to ECG_sig
[y_hat_rls_ecg, W] = RLS(x, r, 15, 0.999);

%plot all the outputs
figure;
subplot(6,1,1)
plot(t, y_i)
title('Desired ECG signal');xlabel('Time(s)');ylabel('Voltage (mV)');
subplot(6,1,2)
plot(t, x)
title('Noisy ECG signal');xlabel('Time(s)');ylabel('Voltage (mV)');
subplot(6,1,3)
plot(t, y_hat_lms_ecg)
title('Filtered ECG Signal using LMS method M=15 , u=0.001');xlabel('Time(s)');ylabel('Voltage (mV)');
subplot(6,1,4)
plot(t, abs(y_hat_lms_ecg - y_i))
title('Absolute Error (LMS)');xlabel('Time(s)');ylabel('Voltage (mV)');
subplot(6,1,5)
plot(t,y_hat_rls_ecg);
title('Filtered ECG Signal using RLS method M=15 lambda = 0.999'); xlabel('Time(s)');ylabel('Voltage (mV)');
subplot(6,1,6)
plot(t,abs(y_hat_rls_ecg - y_i))
title('Absolute Error (RLS)');xlabel('Time(s)');ylabel('Voltage (mV)');
linkaxes();
