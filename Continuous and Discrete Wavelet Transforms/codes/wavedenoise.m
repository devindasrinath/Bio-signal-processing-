% w_name : wavelet name
% y : noisy signal
% num_levels : levels of dcomposition
% th : threshold
function [x_rec ,c_sorted_decend] = wavedenoise( w_name,y, num_levels, th)

    [c, l] = wavedec(y, num_levels, w_name); %decomposition of the 1-D signal y at level num_levels using the wavelet wname(w_name).
    
    c_sorted_decend = sort(abs(c(:)),'descend'); % sorting all the coefficients in decending order in single row vector

    %set coeffcents to zero whch are less than to threhsold
    c_filtered = c;
    c_filtered(abs(c)<th) = 0;

    % reconstruct the signal with the filtered coefficients
    x_rec = waverec(c_filtered, l,w_name );
end