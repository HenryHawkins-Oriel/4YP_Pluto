function [BER_avrg] = BER_vs_Gain(rxSig_BER, avrg_type, sim_params)
    % Characterise Alice and Bob's channel using a plot of BER vs Tx Gain

    % Preallocated vectors for plotting
    vec_template = zeros(1, length(sim_params.tx_gain), length(sim_params.chip_no));
    BER_avrg = vec_template;
    lines_plot = zeros(length(sim_params.tx_gain), length(sim_params.chip_no));
    colours = lines(length(sim_params.chip_no));

    for z = 1 : length(sim_params.chip_no)
    
        % Convert BER cell array into more acessible matrix 
        % Indexes column-by-column
        rxSig_BER_mat = [rxSig_BER{:,:,z}];
        % First row of matrix is BER values
        % Reshape so each column corresponds to one value of tx_gain
        rxSig_BER_mat_shaped = reshape(rxSig_BER_mat(1,:), [], length(sim_params.tx_gain));
        
        % Averaging of SNR values for data repetitions
        % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
        % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
        if strcmp(avrg_type,'mean')
            BER_avrg(1,:,z) = mean(rxSig_BER_mat_shaped, 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            BER_avrg(1,:,z) = median(rxSig_BER_mat_shaped, 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end
        
        % Plot average BER values
        % MATLAB will ignore NaN values in a plot
        figure(1)
        lines_plot(:,z) = plot(sim_params.tx_gain, BER_avrg(1,:,z), ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',15, ...
            'DisplayName', num2str(sim_params.chip_no(z)));
        grid on
        title("BER vs Tx Gain for Spread " + sim_params.mod_type + " at " + num2str(sim_params.f_c/1e6) + ' MHz');
        xlabel('Transmitter Gain (dB)');
        ylabel('Bit Error Rate (s^{-1})');
        ylim([0 10])
        hold on
        
        % Plot scattered (raw) BER values
        tx_gain_vec = repelem(sim_params.tx_gain, 1, sim_params.data_reps);
        scatter(tx_gain_vec, rxSig_BER_mat(1,:), 50, colours(z,:), "x")
        set(gca, 'YScale', 'log')
        hold on;

    end

    hold off
    % Below is temperamental
    for x = 1 : length(sim_params.tx_gain)
        lgd = legend(lines_plot(x,:));
    end
    title(lgd,'Spreading Factor')
    
end