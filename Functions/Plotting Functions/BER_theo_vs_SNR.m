function [SNR_dB_avrg_sort, BER_theo_avrg_sort] = BER_theo_vs_SNR(SNR_dB, BER_theo, avrg_type, sim_params)
    % Characterise Alice and Bob's channel using a plot of BER vs SNR
    
    % Averaging of SNR values for data repetitions
    % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
    % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
    if strcmp(avrg_type,'mean')
        % Mean of each column, across pages
        SNR_dB_avrg = mean(SNR_dB, [1 3],'omitnan');
        BER_theo_avrg = mean(BER_theo, [1 3],'omitnan');
    elseif strcmp(avrg_type,'median')
        SNR_dB_avrg = median(SNR_dB, [1 3],'omitnan');
        BER_theo_avrg = median(BER_theo, [1 3],'omitnan');
    else
        disp("Please input a valid average type.")
    end
    
    % Sort results for line plot
    [SNR_dB_avrg_sort, SNR_index] = sort(SNR_dB_avrg);
    BER_theo_avrg_sort = BER_theo_avrg(SNR_index);
    
    % Plot average BER values
    % MATLAB will ignore NaN values in a plot
    figure(1)
    plot(SNR_dB_avrg_sort, BER_theo_avrg_sort, ...
        'LineStyle','-', ...
        'Marker','.', ...
        'Color', lines(1), ...
        'MarkerSize',15);
    grid on
    title("Theoretical BER vs SNR for Spread " + sim_params.mod_type + " at " + num2str(sim_params.f_c/1e6) + ' MHz');
    xlabel('Signal-to-Noise Ratio (dB)');
    ylabel('Theoretical Bit Error Rate (s^{-1})');
    set(gca, 'YScale', 'log')
    ylim([0 10])
    hold off

    %{
    % Plot scattered (raw) BER values
    % Raw SNR data for scatter plot
    SNR_dB_vec = reshape(SNR_dB,1,[]);
    BER_theo_vec = reshape(BER_theo,1,[]);
    scatter(SNR_dB_vec, BER_theo_vec, 50, lines(1), "x")
    set(gca, 'YScale', 'log')
    hold off
    %}
    
end