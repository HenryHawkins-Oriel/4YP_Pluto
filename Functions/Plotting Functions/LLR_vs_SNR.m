function [A,B,C] = LLR_vs_SNR(SNR_dB, LLR_noise, LLR_sig, results_type, params)
    % Two line plots will allow us to find an accurate threshold in the middle

    % Power of noise shouldn't change
    % Because receiver sensitivity remains the same
    % Therefore, threshold shouldn't change
    % Power of received signal increases with tx_gain
    % Therefore, so does LLR sum (which is proportional to avrg signal power)
    % Which makes it more likely that Willie will detect

    % Choose between simulation or practical results
    if strcmp(results_type,'sim')
        gain_length = length(params.SNR_dB);
    elseif strcmp(results_type,'real')
        gain_length = length(params.tx_gain);
    else
        disp("Please input a valid result type.")
    end

    % Preallocated threshold parameter vector
    A = zeros(2, length(params.chip_no));
    B = zeros(2, length(params.chip_no));
    C = zeros(2, length(params.chip_no));
    % Preallocated vectors for plotting
    lines_plot = zeros(gain_length, length(params.chip_no), 3);
    colours = distinguishable_colors(length(params.chip_no));
    % Declare new figure before looping
    figure
    
    for z = 1 : length(params.chip_no)
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
        % Starting points might need to change with SF...
        % ...and depend on SNR region of interest
        A(:,z) = fminsearch(cost_fcn_a, [0;10], options); % [1;10]
        B(:,z) = fminsearch(cost_fcn_b, [0;10], options); % [1;10]
        C(:,z) = fminsearch(cost_fcn_c, [0;10], options); % [1;10]
        
        % Plot fit functions
        % MATLAB will ignore NaN values in a plot
        % LLR when only noise is present
        % LLR when signal is present amongst noise
        %{
        lines_plot(:,z,1) = plot(SNR_dB_range, fit_fcn(A(:,z), SNR_dB_range), ...
            'LineStyle','-', ...
            'Color', colours(z,:), ...
            'DisplayName', strcat(num2str(params.chip_no(z))));
        hold on
        %lines_plot(:,z,2) = plot(SNR_dB_range, fit_fcn(C(:,z), SNR_dB_range), ...
            'LineStyle','-', ...
            'Color', colours(z,:), ...
            'DisplayName', strcat(num2str(params.chip_no(z))));
        hold on
        %}
        lines_plot(:,z,3) = plot(SNR_dB_range, fit_fcn(B(:,z), SNR_dB_range), ...
            'LineStyle','--', ...
            'Color', colours(z,:), ...
            'DisplayName', strcat(num2str(params.chip_no(z))));
        hold on
        set(gca, 'YScale', 'log')
        grid on
        % Choose title for simulation or practical results
        if strcmp(results_type,'sim')
            %title("LLR vs SNR using Simulated Spread " + params.mod_type);
        elseif strcmp(results_type,'real')
            %title("LLR vs SNR using Spread " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
        else
            disp("Please input a valid result type.")
        end
        xlabel('Signal-to-Noise Ratio (dB)');
        ylabel('Log-Likelihood Ratio');
        xlim([-50,-20]) % high [-15,20] or low [-55,0]
        ylim([1e6,1e12]) % high [1e8,1e13] or low [1e6,1e12]
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
        %scatter(SNR_dB_vec, LLR_sig_noise_avrg, 10, colours(z,:), "o")
        %set(gca, 'YScale', 'log')
        hold on

    end

    hold off
    % Legend
    for x = 1 : gain_length
        lgd = legend(lines_plot(x,:,3));
    end
    title(lgd,'Spreading Factor')
    lgd.Location = 'northeastoutside';
    
end