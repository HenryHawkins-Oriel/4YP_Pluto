function [SNR_dB_avrg_sort, BER_avrg_sort] = BER_vs_SNR(rxSig_BER, SNR_dB, avrg_type, sim_params)
    % Characterise Alice and Bob's channel using a plot of BER vs SNR

    % Preallocated vectors for plotting
    vec_template = zeros(1, length(sim_params.tx_gain), length(sim_params.chip_no));
    SNR_dB_avrg_sort = vec_template;
    BER_avrg_sort = vec_template;
    lines_plot = zeros(length(sim_params.tx_gain), length(sim_params.chip_no));
    colours = lines(length(sim_params.chip_no));

    for z = 1 : length(sim_params.chip_no)
    
        % Convert BER cell array into more acessible matrix 
        % Indexes column-by-column
        rxSig_BER_mat = [rxSig_BER{:,:,z}];
        % First row of matrix is BER values
        % Reshape so each column corresponds to one value of tx_gain
        rxSig_BER_mat_shaped = reshape(rxSig_BER_mat(1,:), sim_params.data_reps, length(sim_params.tx_gain));
        
        % Averaging of SNR values for data repetitions
        % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
        % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
        if strcmp(avrg_type,'mean')
            SNR_dB_avrg = mean(SNR_dB(:,:,z), 1,'omitnan');
            BER_avrg = mean(rxSig_BER_mat_shaped, 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            SNR_dB_avrg = median(SNR_dB(:,:,z), 1,'omitnan');
            BER_avrg = median(rxSig_BER_mat_shaped, 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end
        
        % Sort results for line plot
        [SNR_dB_avrg_sort(1,:,z), SNR_index] = sort(SNR_dB_avrg);
        BER_avrg_sort(1,:,z) = BER_avrg(SNR_index);
        
        % Plot average BER values
        % MATLAB will ignore NaN values in a plot
        figure(1)
        lines_plot(:,z) = plot(SNR_dB_avrg_sort(1,:,z), BER_avrg_sort(1,:,z), ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',15, ...
            'DisplayName', num2str(sim_params.chip_no(z)));
        grid on
        title("BER vs SNR for Spread " + sim_params.mod_type + " at " + num2str(sim_params.f_c/1e6) + ' MHz');
        xlabel('Signal-to-Noise Ratio (dB)');
        ylabel('Bit Error Rate (s^{-1})');
        set(gca, 'YScale', 'log')
        ylim([0 1])
        hold on
        
        %{
        % Plot scattered (raw) BER values
        % Raw SNR data for scatter plot
        SNR_dB_vec = reshape(SNR_dB(:,:,z),1,[]);
        scatter(SNR_dB_vec, rxSig_BER_mat(1,:), 50, colours(z,:), "x")
        set(gca, 'YScale', 'log')
        hold on;
        %}

    end

    hold off
    % Legend
    for x = 1 : length(sim_params.tx_gain)
        lgd = legend(lines_plot(x,:));
    end
    title(lgd,'Spreading Factor')
    
end