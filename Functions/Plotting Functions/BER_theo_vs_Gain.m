function [BER_theo_avrg] = BER_theo_vs_Gain(BER_theo, avrg_type, params)
    % Characterise Alice and Bob's channel using a plot of BER vs Tx Gain

    % Preallocated vectors for plotting
    lines_plot = zeros(length(params.tx_gain), length(params.chip_no));
    colours = distinguishable_colors(length(params.chip_no));

    for z = 1 : length(params.chip_no)

        % We haven't taken spreading into account when calculating BER
        % So we will take it into account here instead
        % Iterate over all elements and pages of BER_theo
        n = params.chip_no(z);
        for j = 1 : length(params.tx_gain)
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

        % Averaging of theoretical BER values
        if strcmp(avrg_type,'mean')
            BER_theo_avrg = mean(BER_theo(:,:,z), 1,'omitnan');
        elseif strcmp(avrg_type,'median')
            BER_theo_avrg = median(BER_theo(:,:,z), 1,'omitnan');
        else
            disp("Please input a valid average type.")
        end

        % Plot average theoretical BER values
        % MATLAB will ignore NaN values in a plot
        figure
        lines_plot(:,z) = plot(params.tx_gain, BER_theo_avrg, ...
            'LineStyle','-', ...
            'Marker','.', ...
            'Color', colours(z,:), ...
            'MarkerSize',6, ...
            'DisplayName', num2str(params.chip_no(z)));
        grid on
        %title("Theoretical BER vs Transmitter Gain using Spread " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
        xlabel('Transmitter Gain (dB)');
        ylabel('Theoretical Bit Error Ratio');
        set(gca, 'YScale', 'log')
        ylim([0 1])
        hold on
   
    
        %{
        % Plot scattered (raw) theoretical BER values
        % Raw data for scatter plot
        tx_gain_mat = repmat(params.tx_gain, params.data_reps, length(params.chip_no));
        tx_gain_vec = reshape(tx_gain_mat,1,[]);
        BER_theo_vec = reshape(BER_theo,1,[]);
        scatter(tx_gain_vec, BER_theo_vec, 50, colours(z,:), "x")
        set(gca, 'YScale', 'log')
        hold off
        %}

    end

    hold off
    % Legend
    for x = 1 : length(params.tx_gain)
        lgd = legend(lines_plot(x,:));
    end
    title(lgd,'Spreading Factor')
    
end