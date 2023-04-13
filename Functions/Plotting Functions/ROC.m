function [TPR, FPR, b] = ROC_v2(A, B, C, LLR_noise, LLR_sig, SNR_dB, sim_params)
    % Plots ROC curve
    % Doesn't change gradient

    % Initialise result vectors
    intcpt_vec_len = 1e3;
    b = zeros(2, intcpt_vec_len);
    TPR = zeros(intcpt_vec_len, length(sim_params.tx_gain), length(sim_params.chip_no));
    FPR = zeros(intcpt_vec_len, length(sim_params.tx_gain), length(sim_params.chip_no));
    H0_sig = zeros(sim_params.data_reps, intcpt_vec_len, length(sim_params.tx_gain), length(sim_params.chip_no));
    H0_noise = zeros(sim_params.data_reps, intcpt_vec_len, length(sim_params.tx_gain), length(sim_params.chip_no));
    % Preallocated vectors for plotting
    lines_plot = zeros(intcpt_vec_len, length(sim_params.tx_gain));
    colours = lines(length(sim_params.tx_gain));

    % For each new spreading sequence length, I get a new figure
    for z = 1 : length(sim_params.chip_no)

        % Curve fit parameters change with number of chips
        a = A(:,z);
        c = C(:,z);
        % Ranges of intercepts
        intcpt_range = c(2) - a(2);
        % Dividing factor dictates range of thresholds
        intcpt_min = a(2) - (intcpt_range/10);
        intcpt_max = c(2) + (intcpt_range/10);
        % Create gradient and intercept vectors for iterating
        grad = B(1);
        intcpt = linspace(intcpt_min, intcpt_max, intcpt_vec_len);

        % For each new tx_gain, I get a new line (or curve)
        for x = 1 : length(sim_params.tx_gain)
            % For each new set (pair) of threshold parameters, ... 
            % I get a new pair of probability values to be plotted
            for i = 1 : length(grad)
                for j = 1 : length(intcpt)
                    idx = (i - 1) * length(intcpt) + j;
                    % Construct threshold parameter vector
                    b(:,idx) = [grad(i); intcpt(j)];
                    % Willie's performance statistics
                    [TPR(idx,x,z), FPR(idx,x,z), H0_sig(:,idx,x,z), H0_noise(:,idx,x,z)] = willie_performance_v3(b(:,idx), LLR_noise(:,x,z), LLR_sig(:,x,z), SNR_dB(:,x,z), sim_params);
                end
            end

            % Sort results for line plot
            [FPR_sort, FPR_index] = sort(FPR(:,x,z));
            TPR_temp = TPR(:,x,z);
            TPR_sort = TPR_temp(FPR_index);

            figure(z)
            lines_plot(:,x) = plot(FPR_sort, TPR_sort, ...
                'LineStyle','-', ...
                'Marker','.', ...
                'Color', colours(x,:), ...
                'MarkerSize',10, ...
                'DisplayName', strcat('Gain: ', num2str(sim_params.tx_gain(x)),' dB'));
            ylim([0 1])
            xlim([0 1])
            grid on
            title("TPR vs FPR with SF " + sim_params.chip_no(z) + " using " + sim_params.mod_type + " at " + num2str(sim_params.f_c/1e6) + ' MHz');
            xlabel('False Positive Rate');
            ylabel('True Positive Rate');
            hold on
            
            % Plot scattered (raw) LLR values
            % LLR when only noise is present
            scatter(FPR(:,x,z), TPR(:,x,z), 10, colours(x,:), "x")
            ylim([0 1])
            xlim([0 1])
            hold on

        end

        hold off
        % Legend
        for idx = 1 : intcpt_vec_len
            lgd = legend(lines_plot(idx,:));
        end
        title(lgd,'Transmitter Gain')

    end

end