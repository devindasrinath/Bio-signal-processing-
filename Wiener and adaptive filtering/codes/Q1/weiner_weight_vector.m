% y : signal desired
% noise : noise of the signal
% M : M of the filter
function W = weiner_weight_vector(y, noise, M)
    % for autocorrelation
    yy_T = 0; 
    nn_T = 0;
    
    % for crosscorrelation
    Yy = 0;
    y_mat = zeros(M,1);
    n_mat = zeros(M,1);

    for i=1:length(y)
        
        y_mat(1) = y(i); 
        n_mat(1) = noise(i);

        %calculate auto correlations
        yy_T = yy_T + toeplitz(autocorr(y_mat, M-1));
        nn_T = nn_T + toeplitz(autocorr(n_mat, M-1));

        %calculate cross correlations
        Yy = Yy + y_mat*y(i);

        % shifting the delay 
        y_mat(2:M) = y_mat(1 : M-1);
        n_mat(2:M) = n_mat(1 : M-1);
    end

    yy_T = yy_T.*mean(y.^2);
    nn_T = nn_T.*mean(noise.^2);

    autocorr_Y = yy_T./ (length(y) - M);
    autocorr_N = nn_T./ (length(y) - M);
    crosscorr_Yy = Yy./ (length(y)-M);
    
    autocorr_X = autocorr_Y + autocorr_N;
    W = autocorr_X\crosscorr_Yy;

end