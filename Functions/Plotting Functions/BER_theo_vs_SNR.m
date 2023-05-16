function [SNR_dB_avrg_sort, BER_theo_avrg_sort] = BER_theo_vs_SNR(SNR_dB, BER_theo, avrg_type, results_type, params)
    % Characterise Alice and Bob's channel using a plot of theoretical BER vs SNR

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
    SNR_dB_avrg_sort = vec_template;
    BER_theo_avrg_sort = vec_template;
    lines_plot = zeros(gain_length, length(params.chip_no));
    colours = distinguishable_colors(length(params.chip_no));
    % Declare new figure before looping
    figure
    
    for z = 1 : length(params.chip_no)
        
        % We haven't taken spreading into account when calculating BER
        % So we will take it into account here instead
        % Iterate over all elements and pages of BER_theo
        n = params.chip_no(z);
        for j = 1 : gain_length
            for i = 1 : params.data_reps
                % Initialise sum of probabilities
                bin_sum = 0;
                if rem(n,2) == 0 % Even SF 
                    for idx = n/2 : n
                        % Binomial coefficient
                        % Number of ways channel can introduce idx errors
                        bin_coeff = nchoosek(n,idx);
                        % Add probability that idx transitions occur...
                        % ...and n−idx bits are received correctly
                        bin_sum = bin_sum + bin_coeff*(BER_theo(i,j,z)^idx)*((1 - BER_theo(i,j,z))^(n - idx));
                    end
                    % Assign value to new BER_theo value
                    BER_theo(i,j,z) = bin_sum;
                else % Odd SF
                    for idx = ((n + 1)/2) : n
                        % Binomial coefficient
                        % Number of ways channel can introduce idx errors
                        bin_coeff = nchoosek(n,idx);
                        % Add probability that idx transitions occur...
                        % ...and n−idx bits are received correctly
                        bin_sum = bin_sum + bin_coeff*(BER_theo(i,j,z)^idx)*((1 - BER_theo(i,j,z))^(n - idx));
                    end
                    % Assign value to new BER_theo value
                    BER_theo(i,j,z) = bin_sum;
                end
            end
        end
    
        % Averaging of SNR values for data repetitions
        % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
        % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
        if strcmp(avrg_type,'mean')
            SNR_dB_avrg = mean(SNR_dB(:,:,z), 1,'omitnan');
            BER_theo_avrg = mean(BER_theo(:,:,z), 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            SNR_dB_avrg = median(SNR_dB(:,:,z), 1,'omitnan');
            BER_theo_avrg = median(BER_theo(:,:,z), 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end
        
        % Sort results for line plot
        [SNR_dB_avrg_sort(1,:,z), SNR_index] = sort(SNR_dB_avrg);
        BER_theo_avrg_sort(1,:,z) = BER_theo_avrg(SNR_index);

        % Plot average BER values
        % MATLAB will ignore NaN values in a plot
        lines_plot(:,z) = plot(SNR_dB_avrg_sort(1,:,z), BER_theo_avrg_sort(1,:,z), ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',6, ...
            'DisplayName', num2str(params.chip_no(z)));
        grid on
        % Choose title for simulation or practical results
        if strcmp(results_type,'sim')
            %title("Theoretical BER vs SNR using Simulated Spread " + params.mod_type);
        elseif strcmp(results_type,'real')
            %title("Theoretical BER vs SNR using Spread " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
        else
            disp("Please input a valid result type.")
        end
        xlabel('Signal-to-Noise Ratio (dB)');
        ylabel('Theoretical Bit Error Ratio');
        set(gca, 'YScale', 'log')
        ylim([1e-4,1])
        xlim([-30,15])
        hold on
   
    
        %{
        % Plot scattered (raw) theoretical BER values
        % Raw SNR data for scatter plot
        SNR_dB_vec = reshape(SNR_dB,1,[]);
        BER_theo_vec = reshape(BER_theo,1,[]);
        scatter(SNR_dB_vec, BER_theo_vec, 50, colours(z,:), "x")
        set(gca, 'YScale', 'log')
        hold off
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