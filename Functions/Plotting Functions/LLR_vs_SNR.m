function [A,B,C] = LLR_vs_SNR(SNR_dB, LLR_noise, LLR_sig, avrg_type, sim_params)
    % Two line plots will allow us to find an accurate threshold in the middle

    % Power of noise shouldn't change
    % Because receiver sensitivity remains the same
    % Therefore, threshold shouldn't change
    % Power of received signal increases with tx_gain
    % Therefore, so does LLR sum (which is proportional to avrg signal power)
    % Which makes it more likely that Willie will detect

    % Preallocated threshold parameter vector
    A = zeros(2, length(sim_params.chip_no));
    B = zeros(2, length(sim_params.chip_no));
    C = zeros(2, length(sim_params.chip_no));
    % Preallocated vectors for plotting
    lines_plot = zeros(length(sim_params.tx_gain), length(sim_params.chip_no), 3);
    colours = lines(length(sim_params.chip_no));

    for z = 1 : length(sim_params.chip_no)
        %{
        % Averaging of SNR values for data repetitions
        % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
        % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
        if strcmp(avrg_type,'mean')
            SNR_dB_avrg = mean(SNR_dB(:,:,z), 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            SNR_dB_avrg = median(SNR_dB(:,:,z), 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end
        
        % Sort results for line plot
        [SNR_dB_avrg_sort, SNR_index] = sort(SNR_dB_avrg);
        
        % Averaging of LLR values from data repetitions
        if strcmp(avrg_type,'mean')
            LLR_noise_avrg = mean(LLR_noise(:,:,z), 1,'omitnan');
            LLR_sig_avrg = mean(LLR_sig(:,:,z), 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            LLR_noise_avrg = median(LLR_noise(:,:,z), 1,'omitnan');
            LLR_sig_avrg = median(LLR_sig(:,:,z), 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end
        

        % Sort results for line plot
        LLR_noise_avrg_sort = LLR_noise_avrg(SNR_index);
        LLR_sig_avrg_sort = LLR_sig_avrg(SNR_index);
        %}

        % Raw SNR data for scatter plot (remove NaNs)
        SNR_dB_vec = reshape(SNR_dB(:,:,z),1,[]);
        SNR_dB_vec_no_NaNs = SNR_dB_vec(~isnan(SNR_dB_vec));
        % Raw LLR data for scatter plot (remove NaNs)
        LLR_noise_vec = reshape(LLR_noise(:,:,z),1,[]);
        LLR_noise_vec_no_NaNs = LLR_noise_vec(~isnan(LLR_noise_vec));
        LLR_sig_vec = reshape(LLR_sig(:,:,z),1,[]);
        LLR_sig_vec_no_NaNs = LLR_sig_vec(~isnan(LLR_sig_vec));
        % Decision boundary taken to be values exactly halfway between the above
        LLR_sig_noise_avrg = LLR_noise_vec + 0.5*(LLR_sig_vec - LLR_noise_vec);
        LLR_sig_noise_avrg_no_NaNs = LLR_sig_noise_avrg(~isnan(LLR_sig_noise_avrg));
        
        % Range over which to plot fit functions
        SNR_dB_range = linspace(min(SNR_dB_vec), max(SNR_dB_vec));
    
        % Curve fitting
        % Objective Function
        fit_fcn = @(t,x) 10.^(t(1)*x + t(2));
        % Cost functions are 2-norms of differences
        % norm() does not accept NaN values as inputs
        cost_fcn_a = @(a) norm(LLR_sig_vec_no_NaNs - fit_fcn(a,SNR_dB_vec_no_NaNs));
        cost_fcn_b = @(b) norm(LLR_sig_noise_avrg_no_NaNs - fit_fcn(b,SNR_dB_vec_no_NaNs));
        cost_fcn_c = @(c) norm(LLR_noise_vec_no_NaNs - fit_fcn(c,SNR_dB_vec_no_NaNs));
        % Estimate parameters
        % Start at (1,10) and find local minimum (A,B,C) of cost functions
        options = optimset('fminsearch');
        options.MaxIter = 1000;
        options.MaxFunEvals = 1000;        
        % Starting values might need to change with z (SF)?
        A(:,z) = fminsearch(cost_fcn_a, [1;10], options);
        B(:,z) = fminsearch(cost_fcn_b, [1;10], options);
        C(:,z) = fminsearch(cost_fcn_c, [1;10], options);
        
        % Plot fit functions
        % MATLAB will ignore NaN values in a plot
        % LLR when only noise is present
        % LLR when signal is present amongst noise
        figure(2)
        lines_plot(:,z,1) = plot(SNR_dB_range, fit_fcn(A(:,z), SNR_dB_range), ...
            'LineStyle','-', ...
            'Color', colours(z,:), ...
            'DisplayName', strcat('Threshold (', num2str(sim_params.chip_no(z)),')'));
        hold on
        lines_plot(:,z,2) = plot(SNR_dB_range, fit_fcn(C(:,z), SNR_dB_range), ...
            'LineStyle','-', ...
            'Color', colours(z,:), ...
            'DisplayName', strcat('Threshold (', num2str(sim_params.chip_no(z)),')'));
        hold on
        lines_plot(:,z,3) = plot(SNR_dB_range, fit_fcn(B(:,z), SNR_dB_range), ...
            'LineStyle','--', ...
            'Color', colours(z,:), ...
            'DisplayName', strcat('Threshold (', num2str(sim_params.chip_no(z)),')'));
        hold on
        set(gca, 'YScale', 'log')
        grid on
        title("LLR vs SNR for Spread " + sim_params.mod_type + " at " + num2str(sim_params.f_c/1e6) + ' MHz');
        %legend()
        xlabel('Signal-to-Noise Ratio (dB)');
        ylabel('Log-Likelihood Ratio (LLR)');
        hold on
        
        % Plot scattered (raw) LLR values
        % LLR when only noise is present
        scatter(SNR_dB_vec, LLR_sig_vec, 10, colours(z,:), "v")
        set(gca, 'YScale', 'log')
        hold on
        % LLR when signal is present amongst noise
        scatter(SNR_dB_vec, LLR_noise_vec, 10, colours(z,:), "^")
        set(gca, 'YScale', 'log')
        % LLR threshold (function that fits average of both LLRs)
        scatter(SNR_dB_vec, LLR_sig_noise_avrg, 10, colours(z,:), "o")
        set(gca, 'YScale', 'log')
        hold on

    end

    hold off
    % Legend
    for x = 1 : length(sim_params.tx_gain)
        for y = 1 : 3
            lgd = legend(lines_plot(x,:,y));
        end
    end
    title(lgd,'Spreading Factor')
    
end