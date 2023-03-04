function [y_hat, W_f] = weiner_fil_freq(y_i, noise, x)

    S_yy = abs(fft(y_i,length(x)*2-1)).^2;  % Power of desired signal
    S_NN = abs(fft(noise,length(x)*2-1)).^2;% Power of noise
    S_Xf  = fft(x,length(x)*2-1);           % FT of the noisy singnal
    W_f = S_yy./(S_yy + S_NN);          % Freqescies of  Wiener filter 
    S_Y_hat = W_f.*S_Xf;                % Signal estimate using observation and using Wiener filter 
    y_hat_time = (ifft(S_Y_hat));       % time domain representation of signal
    y_hat = y_hat_time(1: length(x));   % positive range extracted

end