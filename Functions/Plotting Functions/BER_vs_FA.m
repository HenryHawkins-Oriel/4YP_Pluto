function [prob_fa_vec, H0_noise] = BER_vs_FA(b, SNR_dB, rxSig_BER, LLR_noise, avrg_type, results_type, params)
    % Characterise Alice and Willie's channel using a plot of BER vs FA

    % H0 = false means "Transmission detected." (H0 rejected)
    % H0 = true means "Noise detected." (H0 accepted)

    % Choose between simulation or practical results
    if strcmp(results_type,'sim')
        gain_length = length(params.SNR_dB);
    elseif strcmp(results_type,'real')
        gain_length = length(params.tx_gain);
    else
        disp("Please input a valid result type.")
    end

    % Initialise H0 for the noise case
    H0_noise = false(params.data_reps, gain_length, length(params.chip_no));
    H0_noise_grnd_truth = true(params.data_reps, gain_length, length(params.chip_no));

    % If an SNR value is NaN, then we have an invalid test
    for x = 1 : gain_length
        for y = 1 : params.data_reps
            for z = 1 : length(params.chip_no)
                if ~isnan(SNR_dB(y,x,z))
                    % Compute threshold using experimental parameters
                    thresh = 10^(b(1,z)*SNR_dB(y,x,z) + b(2,z));
                    H0_noise(y,x,z) = willie_decision(LLR_noise(y,x,z), thresh);
                else
                    H0_noise(y,x,z) = false;
                end
            end
        end
    end
    
    % NaN value in LLR matrix means an invalid test and a false in H0 matrix
    % Which means H0 was wrongly rejected (false alarm)
    % This false value will contribute to the fa_cnt
    % This false value will also be included in the number of valid tests
    % Subtract number of NaNs from fa_cnt and total number of tests
    LLR_noise_NaNs = sum(isnan(LLR_noise), 1);

    % FALSE POSITIVE RATE (FPR)
    % Calculate probability of false alarm
    % Wrongly rejected H0 when H0 was true
    % Ground truth is matrix of "trues" (accepting H0)
    % Therefore, interested when there is a "false" value
    fa_mat = H0_noise_grnd_truth - H0_noise;
    fa_cnt_vec = abs(sum(fa_mat, 1)) - LLR_noise_NaNs;
    valid_tests_vec = repelem(params.data_reps, 1, gain_length, length(params.chip_no));
    valid_noise_tests = (valid_tests_vec - LLR_noise_NaNs);
    prob_fa_vec = fa_cnt_vec ./ valid_noise_tests;

    % Preallocated plotting vectors
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
        rxSig_BER_mat_shaped = reshape(rxSig_BER_mat(1,:), params.data_reps, gain_length);
        
        % Averaging of SNR values for data repetitions
        if strcmp(avrg_type,'mean')
            BER_avrg = mean(rxSig_BER_mat_shaped, 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            BER_avrg = median(rxSig_BER_mat_shaped, 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end
        
        % Sort results for line plot
        [prob_fa_sort, fa_index] = sort(prob_fa_vec(1,:,z));
        BER_avrg_sort = BER_avrg(fa_index);
        
        % Plot average BER values
        % MATLAB will ignore NaN values in a plot
        lines_plot(:,z) = plot(prob_fa_sort, BER_avrg_sort, ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',15, ...
            'DisplayName', num2str(params.chip_no(z)));
        grid on
        % Choose title for simulation or practical results
        if strcmp(results_type,'sim')
            %title("BER vs FPR using Simulated Spread " + params.mod_type);
        elseif strcmp(results_type,'real')
            %title("BER vs FPR using Spread " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
        else
            disp("Please input a valid result type.")
        end
        xlabel('False Positive Rate');
        ylabel('Bit Error Ratio');
        set(gca, 'YScale', 'log')
        xlim([0,1])
        ylim([1e-4,1])
        hold on
    
    end
    
    hold off
    % Legend
    for x = 1 : gain_length
        lgd = legend(lines_plot(x,:));
    end
    title(lgd,'Spreading Factor')
    lgd.Location = 'southeast';
    
end