%%INPUTS
% x : noisy signal (row vector)
% r : referance signal (row vector)
% M : filter length
% mu : rate of convergnace

%%OUTPUTS
% y_hat : filtered signal
% W : final weigths
function [y_hat, W] = LMS(x , r, M, mu)

    % exception handling
    % for array sizes
    if (length(x)<= M)  
        disp('error: Signal length is less than the filter order, please reduce the order');
        return; 
    end
    if (length(x) ~= length(r))  
        disp('error: Input signal and reference signal are not the same size');
        return; 
    end
    
    %for mu
    lambda_max = 20*M*((x*x')/length(x));
    if mu > 2/lambda_max 
        disp(['mu is too large' num2str(mu) ' >' num2str(2/lambda_max) , ' please reduce mu.']);
        return
    end

    W = zeros(M ,1); % W have M number of rows (w(0) ; w(1).....; w(M-1))
    r_padded = [zeros(1,M-1) r]; %padding M-1 number of zeros at the start
    
    e = 0; % estimate(error of ANC)
    
    for i = 1 : length(x)
        R = (flip(r_padded(i:M+i-1)))';% R have M number of rows (R(M-1) ; R(M-2).....; R(0))
        e = x(i) - W'*R; %calculate the error
        W = W +2*mu*e*R;
        y_hat(i) = e; %system output
    end
    
    

end