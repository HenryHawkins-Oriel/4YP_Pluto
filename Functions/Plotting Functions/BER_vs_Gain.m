function [BER_avrg] = BER_vs_Gain(rxSig_BER, avrg_type, params)
    % Characterise Alice and Bob's channel using a plot of BER vs Tx Gain

    % Preallocated vectors for plotting
    vec_template = zeros(1, length(params.tx_gain), length(params.chip_no));
    BER_avrg = vec_template;
    lines_plot = zeros(length(params.tx_gain), length(params.chip_no));
    colours = distinguishable_colors(length(params.chip_no));
    % Declare new figure before looping
    figure
    
    for z = 1 : length(params.chip_no)
    
        % Convert BER cell array into more acessible matrix 
        % Indexes column-by-column
        rxSig_BER_mat = [rxSig_BER{:,:,z}];
        % First row of matrix is BER values
        % Reshape so each column corresponds to one value of tx_gain
        rxSig_BER_mat_shaped = reshape(rxSig_BER_mat(1,:), [], length(params.tx_gain));
        
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
        figure(3)
        lines_plot(:,z) = plot(params.tx_gain, BER_avrg(1,:,z), ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',15, ...
            'DisplayName', num2str(params.chip_no(z)));
        grid on
        %title("BER vs Transmitter Gain for Spread " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
        xlabel('Transmitter Gain (dB)');
        ylabel('Bit Error Ratio');
        ylim([1e-4 1])
        hold on
        
        % Plot scattered (raw) BER values
        tx_gain_vec = repelem(params.tx_gain, 1, params.data_reps);
        scatter(tx_gain_vec, rxSig_BER_mat(1,:), 50, colours(z,:), "x")
        set(gca, 'YScale', 'log')
        hold on;

    end

    hold off
    % Below is temperamental
    for x = 1 : length(params.tx_gain)
        lgd = legend(lines_plot(x,:));
    end
    title(lgd,'Spreading Factor')
    
end