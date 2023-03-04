% w_name : wavlet name
% x_name : signal name of expected
% y_name : signal name of expected
% coeff : wavelt coeffcents(approx. , detailed)
% x_rec : reconstructed signal
% x : orginal signal
function plot_denoise_wave(w_name ,x_name ,y_name, coeff, x_rec, x)

    % Plotthe sorted coefficients
    figure;
    stem(coeff);
    title(['Sorted ' w_name ' wavelet coefficients of  ' y_name ' | Descending order']);
    xlim([0, length(coeff)]);
    
    %plot reconstruted signal
    figure;
    plot(1:length(x_rec), x_rec);
    title([x_name ' reconstructed  |' w_name]), xlabel('Samples(n)'), ylabel('Amplitude');
    xlim([0 length(x_rec)]);
    
    % Compare original and the reconstructed signal using plots
    figure;
    plot(1 : length(x_rec), x, 1 : length(x_rec), x_rec)
    title(['Comparing original and reconstructed ' x_name ' | ' w_name ]), xlabel('Samples(n)'), ylabel('Amplitude');
    legend(x_name, ['Reconstructed ' x_name]);
    xlim([0 length(x_rec)]);
    
end