function [prob_md_vec, H0_sig] = MD_vs_SNR(b, SNR_dB, LLR_sig, avrg_type, results_type, params)
    % Characterise Alice and Willie's channel using a plot of MD vs SNR

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
    H0_sig = true(params.data_reps, gain_length, length(params.chip_no));
    H0_sig_gt = false(params.data_reps, gain_length, length(params.chip_no));

    % If an SNR value is NaN, then we have an invalid test
    for x = 1 : gain_length
        for y = 1 : params.data_reps
            for z = 1 : length(params.chip_no)
                if ~isnan(SNR_dB(y,x,z))
                    % Compute threshold using experimental parameters
                    thresh = 10^(b(1,z)*SNR_dB(y,x,z) + b(2,z));
                    H0_sig(y,x,z) = willie_decision(LLR_sig(y,x,z), thresh);
                else
                    H0_sig(y,x,z) = true;
                end
            end
        end
    end
    
    % NaN value in SNR matrix means an invalid test and a true in H0 matrix
    % Which means H0 was wrongly accepted (missed detection)
    % This true value will contribute to the md_cnt
    % This true value will also be included in the number of valid tests
    % Subtract number of NaNs from md_cnt and total number of tests
    LLR_sig_NaNs = sum(isnan(LLR_sig), 1);

    % FALSE POSITIVE RATE (FPR)
    % Calculate probability of false alarm
    % Wrongly rejected H0 when H0 was true
    % Ground truth is matrix of "trues" (accepting H0)
    % Therefore, interested when there is a "false" value
    md_mat = H0_sig_gt - H0_sig;
    md_cnt_vec = abs(sum(md_mat, 1)) - LLR_sig_NaNs;
    valid_tests_vec = repelem(params.data_reps, 1, gain_length, length(params.chip_no));
    valid_sig_tests = (valid_tests_vec - LLR_sig_NaNs);
    prob_md_vec = md_cnt_vec ./ valid_sig_tests;

    % Averaging of SNR values for data repetitions
    % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
    % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
    if strcmp(avrg_type,'mean')
        SNR_dB_avrg = mean(SNR_dB, 1,'omitnan');
    elseif strcmp(avrg_type,'median')
        SNR_dB_avrg = median(SNR_dB, 1,'omitnan');
    else
        disp("Please input a valid average type.")
    end

    % Preallocated plotting vectors
    lines_plot = zeros(gain_length, length(params.chip_no));
    colours = distinguishable_colors(length(params.chip_no));
    % Declare new figure before looping
    figure

    for z = 1 : length(params.chip_no)
    
        % Sort results for line plot
        [SNR_dB_avrg_sort, SNR_index] = sort(SNR_dB_avrg(1,:,z));
        prob_md_vec_temp = prob_md_vec(1,:,z);
        prob_md_sort = prob_md_vec_temp(SNR_index);
        
        % Plot average BER values
        % MATLAB will ignore NaN values in a plot
        lines_plot(:,z) = plot(SNR_dB_avrg_sort, prob_md_sort, ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',15, ...
            'DisplayName', num2str(params.chip_no(z)));
        grid on
        % Choose title for simulation or practical results
        if strcmp(results_type,'sim')
            %title("FNR vs SNR using Simulated Spread " + params.mod_type);
        elseif strcmp(results_type,'real')
            %title("FNR vs SNR using Spread " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
        else
            disp("Please input a valid result type.")
        end
        legend()
        xlabel('Signal-to-Noise Ratio (dB)');
        ylabel('False Negative Rate');
        ylim([0 1])
        xlim([-15,20])
        hold on

    end

    hold off
    % Legend
    for x = 1 : gain_length
        lgd = legend(lines_plot(x,:));
    end
    title(lgd,'Spreading Factor')
    lgd.Location = 'northeast';
    
end