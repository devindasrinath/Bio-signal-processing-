% w_name : wavelet name
% x : signal
% num_levels : levels of dcomposition
% percentage : percentage of eenrgy to preserve

function [x_compressed ,c_sorted_decend] = wavecompress(w_name, x, num_levels, percentage)

    [c, l] = wavedec(x, num_levels, w_name); %decomposition of the 1-D signal y at level num_levels using the wavelet wname(w_name).
    
    c_sorted_decend = sort(abs(c(:)),'descend'); % sorting all the coefficients in decending order in single row vector
    
    

    E_x = c_sorted_decend'*c_sorted_decend; %total energy of wave
    EC_x = 0; %cummilative energy 
    num_coef_req = 0; %number of coeffcents for 
    for i=1:length(c_sorted_decend)
        EC_x = EC_x + (c_sorted_decend(i)*c_sorted_decend(i));
        if ((EC_x/E_x) >= percentage) %uncomment this if you want exact 99% energy
        %if (round((EC_x/E_x),2) >= percentage)
            num_coef_req = i;
            disp(['Number coefficients required for 99% of the energy of the signal = ' num2str(num_coef_req)]);
            break;
        end
    end

    %compression ration calculation
    compression_ratio = length(x)/num_coef_req;
    disp(['Compression Ratio = ' num2str(compression_ratio) ' : 1']);

    %set 0 to unuse coeffcents
    c_selected = c;
    c_selected(abs(c)<c_sorted_decend(num_coef_req)) = 0;
  
    % reconstruct the signal using remaining coefficients
    x_compressed = waverec(c_selected, l, w_name);

end