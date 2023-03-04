%%INPUTS
% x : noisy signal (row vector)
% r : referance signal (row vector)
% M : filter length
% lambda : forrgetting factor

%%OUTPUTS
% y_hat : filtered signal
% W : final weigths
function [y_hat, W] = RLS(x , r,  M , lambda)

    W = zeros(M ,1); % W have M number of rows (w(0) ; w(1).....; w(M-1))
    r_padded = [zeros(1,M-1) r]; %padding M-1 number of zeros at the start
    
    P  = 0.01*eye(M); %P(0)

    for i = 1 : length(x)
        %calculate error
        R = (flip(r_padded(i:M+i-1)))';% R have M number of rows (R(M-1) ; R(M-2).....; R(0))
        e = x(i) - W'*R; %calculate the error
        
        %calculate K(n)
        K = (P*R)./(lambda + R'*P*R);
        
        %calculate P(n)
        P = (P - K*R'*P)./lambda;
        
        %calculate W(n)
        W = W  + K*e;
        
        y_hat(i) = e;
    end
   
end