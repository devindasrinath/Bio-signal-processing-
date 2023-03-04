%x : signal without noise
%w : weight vector
function y_hat = weiner_filter(x, W)
%     M = length(W);
%     x_padded= [ x zeros(1,M-1)];
%     for i = 1: length(x)
%         x_extracted = (flip(x_padded(i:M+i-1)));% R have M number of rows (R(M-1) ; R(M-2).....; R(0))
%         y_hat(i) = x_extracted * flip(W);
%     end
    order = length(W);
    y_hat = x;
    signal = x;
    weight_mat = W;
    for i = 1: length(signal) - order
        y_hat(i) = signal(i : i + order - 1) * weight_mat;
    end
end