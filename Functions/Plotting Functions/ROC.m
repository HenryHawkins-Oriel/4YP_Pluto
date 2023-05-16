function [TPR, FPR] = ROC(A, B, C, LLR_noise, LLR_sig, SNR_dB, results_type, params)
    % Plots ROC curve
    % Doesn't change gradient

    % No of times the intercept varies
    intcpt_vec_len = 1e4;

    % Choose between simulation or practical results
    if strcmp(results_type,'sim')
        gain_vec = params.SNR_dB;
        H0_sig = zeros(params.data_reps, intcpt_vec_len, length(gain_vec), length(params.chip_no));
        H0_noise = zeros(params.data_reps, intcpt_vec_len, length(gain_vec), length(params.chip_no));
    elseif strcmp(results_type,'real')
        % If results have been obtained practically, SNR values are...
        % ...erratic so will need putting into bins
        % Specify number of bins to sort SNR values into
        no_bins = length(params.tx_gain);
        % Define matrices to be repopulated with binned values
        % It's possible all elements could end up in one column, ...
        % ...so each column needs to be big enough for this possibility
        fill_mat_1 = NaN(numel(SNR_dB), no_bins, length(params.chip_no));
        fill_mat_2 = NaN(numel(SNR_dB), no_bins, length(params.chip_no));
        fill_mat_3 = NaN(numel(SNR_dB), no_bins, length(params.chip_no));
        % Find which bin each SNR_dB value should go into
        for z =  1 : length(params.chip_no)
            SNR_dB_temp = SNR_dB(:,:,z);
            LLR_sig_temp = LLR_sig(:,:,z);
            LLR_noise_temp = LLR_noise(:,:,z);
            [idx, ~] = discretize(SNR_dB_temp, no_bins);
            %centres = edges(1:no_bins) + 0.5*(edges(2) - edges(1));
            for x = 1 : no_bins
                bin_elems = length(SNR_dB(ismember(idx, x)));
                fill_mat_1(1:bin_elems,x,z) = SNR_dB_temp(ismember(idx, x));
                fill_mat_2(1:bin_elems,x,z) = LLR_sig_temp(ismember(idx, x));
                fill_mat_3(1:bin_elems,x,z) = LLR_noise_temp(ismember(idx, x));
            end
        end
        % Assign filled matrices to new SNR and LLR matrices
        SNR_dB = fill_mat_1;
        LLR_sig = fill_mat_2;
        LLR_noise = fill_mat_3;
        % Average of SNR_dB columns for legend
        % If all data points from a tx_gain value (i.e. a column in SNR_dB) are NaN...
        % ... then mean/medians() will also return NaN (even with 'omitnan' flag)
        SNR_dB_avrg = mean(SNR_dB, 1,'omitnan');
        gain_vec = SNR_dB_avrg;
        % Initialise H0 matrices
        % Number of rows has now changed to accomodate grouping into bins
        [data_reps, ~] = size(SNR_dB);
        H0_sig = zeros(data_reps, intcpt_vec_len, length(gain_vec), length(params.chip_no));
        H0_noise = zeros(data_reps, intcpt_vec_len, length(gain_vec), length(params.chip_no));
    else
        disp("Please input a valid result type.")
    end

    % Initialise result vectors
    b = zeros(2, intcpt_vec_len);
    TPR = zeros(intcpt_vec_len, length(gain_vec), length(params.chip_no));
    FPR = zeros(intcpt_vec_len, length(gain_vec), length(params.chip_no));
    % Preallocated vectors for plotting
    lines_plot = zeros(intcpt_vec_len, length(gain_vec));
    colours = distinguishable_colors(length(gain_vec));

    % For each new spreading sequence length, I get a new figure
    for z = 1 : length(params.chip_no)

        % Curve fit parameters change with number of chips
        a = A(:,z);
        c = C(:,z);
        % Ranges of intercepts
        intcpt_range = c(2) - a(2);
        % Dividing factor dictates range of thresholds
        intcpt_min = a(2) - (intcpt_range/2);
        intcpt_max = c(2) + (intcpt_range/2);
        % Create gradient and intercept vectors for iterating
        grad = B(1);
        intcpt = linspace(intcpt_min, intcpt_max, intcpt_vec_len);
        % For each new tx_gain, I get a new line (or curve)
        for x = 1 : length(gain_vec)
            % For each new set (pair) of threshold parameters, ... 
            % I get a new pair of probability values to be plotted
            for i = 1 : length(grad)
                for j = 1 : length(intcpt)
                    idx = (i - 1) * length(intcpt) + j;
                    % Construct threshold parameter vector
                    b(:,idx) = [grad(i); intcpt(j)];
                    % Willie's operating characteristics
                    [TPR(idx,x,z), FPR(idx,x,z), H0_sig(:,idx,x,z), H0_noise(:,idx,x,z)] = ...
                        MD_FA_calc(b(:,idx), LLR_noise(:,x,z), LLR_sig(:,x,z), SNR_dB(:,x,z));
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
                'MarkerSize',5, ...
                'DisplayName', strcat(num2str(round(gain_vec(z,x)))));
            ylim([0 1])
            xlim([0 1])
            grid on
            % Choose title for simulation or practical results
            if strcmp(results_type,'sim')
                %title("TPR vs FPR with SF " + params.chip_no(z) + " using Simulated " + params.mod_type);
            elseif strcmp(results_type,'real')
                %title("TPR vs FPR with SF " + params.chip_no(z) + " using " + params.mod_type + " at " + num2str(params.f_c/1e6) + ' MHz');
            else
                disp("Please input a valid result type.")
            end
            xlabel('False Positive Rate');
            ylabel('True Positive Rate');
            hold on

        end

        hold off
        % Legend
        for idx = 1 : intcpt_vec_len
            lgd = legend(lines_plot(idx,:));
        end
        title(lgd,'SNR (dB)')
        lgd.Location = 'southeast';

    end

end