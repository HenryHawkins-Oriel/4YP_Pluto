function [BER_theo_avrg] = BER_theo_vs_Gain(BER_theo, avrg_type, sim_params)
    % Characterise Alice and Bob's channel using a plot of BER vs Tx Gain
            
    % Averaging of SNR values for data repetitions
    % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
    % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
    if strcmp(avrg_type,'mean')
        BER_theo_avrg = mean(BER_theo, [1 3],'omitnan');
    elseif strcmp(avrg_type,'median')
        BER_theo_avrg = median(BER_theo, [1 3],'omitnan');
    else
        disp("Please input a valid average type.")
    end
    
    % Plot average BER values
    % MATLAB will ignore NaN values in a plot
    figure(1)
    plot(sim_params.tx_gain, BER_theo_avrg, ...
        'LineStyle','-', ...
        'Marker','.', ...
        'Color', lines(1), ...
        'MarkerSize',15);
    grid on
    title("Theoretical BER vs Tx Gain for Spread " + sim_params.mod_type + " at " + num2str(sim_params.f_c/1e6) + ' MHz');
    xlabel('Transmitter Gain (dB)');
    ylabel('Bit Error Rate (s^{-1})');
    ylim([0 10])
    hold on
    
    % Plot scattered (raw) BER values
    BER_theo_vec = reshape(BER_theo,1,[]);
    tx_gain_mat = repmat(sim_params.tx_gain, sim_params.data_reps, length(sim_params.chip_no));
    tx_gain_vec = reshape(tx_gain_mat,1,[]);
    scatter(tx_gain_vec, BER_theo_vec, 50, lines(1), "x")
    set(gca, 'YScale', 'log')
    hold off

end