function [SNR_dB_avrg_sort, BER_avrg_sort, BER_mat] = BER_vs_SNR(rxSig_BER, SNR_dB, avrg_type, results_type, params)
    % Characterise Alice and Bob's channel using a plot of BER vs SNR

    % Choose between simulation or practical results
    if strcmp(results_type,'sim')
        gain_length = length(params.SNR_dB);
    elseif strcmp(results_type,'real')
        gain_length = length(params.tx_gain);
    else
        disp("Please input a valid result type.")
    end

    % Preallocated vectors for plotting
    vec_template = zeros(1, gain_length, length(params.chip_no));
    BER_mat = NaN(params.data_reps, gain_length, length(params.chip_no));
    SNR_dB_avrg_sort = vec_template;
    BER_avrg_sort = vec_template;
    lines_plot = zeros(gain_length, length(params.chip_no));
    colours = distinguishable_colors(length(params.chip_no));
    % Declare new figure before looping
    figure
    
    for z = 1 : length(params.chip_no)
    
        % Convert BER cell array into more acessible matrix 
        % Indexes column-by-column
        rxSig_BER_mat = [rxSig_BER{:,:,z}];
        % First row of matrix is BER values
        % Reshape so each column corresponds to one value of tx_gain
        BER_mat(:,:,z) = reshape(rxSig_BER_mat(1,:), params.data_reps, gain_length);
        
        % Averaging of SNR values for data repetitions
        % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
        % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
        if strcmp(avrg_type,'mean')
            SNR_dB_avrg = mean(SNR_dB(:,:,z), 1,'omitnan');
            BER_avrg = mean(BER_mat(:,:,z), 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            SNR_dB_avrg = median(SNR_dB(:,:,z), 1,'omitnan');
            BER_avrg = median(BER_mat(:,:,z), 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end
        
        % Sort results for line plot
        [SNR_dB_avrg_sort(1,:,z), SNR_index] = sort(SNR_dB_avrg);
        BER_avrg_sort(1,:,z) = BER_avrg(SNR_index);
        
        % Plot average BER values
        % MATLAB will ignore NaN values in a plot
        lines_plot(:,z) = plot(SNR_dB_avrg_sort(1,:,z), BER_avrg_sort(1,:,z), ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',6, ...
            'DisplayName', num2str(params.chip_no(z)));
        grid on
        % Choose title for simulation or practical results
        if strcmp(results_type,'sim')
            %title("BER vs SNR using Simulated Spread " + params.mod_type);
        elseif strcmp(results_type,'real')
            %title("BER vs SNR using Simulated Spread " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
        else
            disp("Please input a valid result type.")
        end
        xlabel('Signal-to-Noise Ratio (dB)');
        ylabel('Bit Error Ratio');
        set(gca, 'YScale', 'log')
        ylim([1e-4, 1])
        xlim([-10,20])
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
    for x = 1 : gain_length
        lgd = legend(lines_plot(x,:));
    end
    title(lgd,'Spreading Factor')
    lgd.Location = 'southwest';
    
end