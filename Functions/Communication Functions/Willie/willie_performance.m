function [TPR, FPR, H0_sig, H0_noise] = willie_performance(b, LLR_noise, LLR_sig, SNR_dB, sim_params)
    % Rejects or accepts null hypothesis based on threshold input (b)
    % Calculates probability of false alarm and missed detection

    %H0 = false means "Transmission detected." (H0 rejected)
    %H0 = true means "Noise detected." (H0 accepted)

    H0_sig = true(sim_params.data_reps,1);
    H0_noise = false(sim_params.data_reps,1);

    % If an LLR value is NaN, then we have an invalid test
    for y = 1 : sim_params.data_reps
        if ~isnan(SNR_dB(y,1))
            % Compute threshold using experimental parameters
            thresh = 10^(b(1)*SNR_dB(y,1) + b(2));
            H0_sig(y,1) = willie_decision(LLR_sig(y,1), thresh);
            H0_noise(y,1) = willie_decision(LLR_noise(y,1), thresh);
        else
            H0_sig(y,1) = true;
            H0_noise(y,1) = false;
        end
    end

    % Ground truth column vectors
    H0_noise_grnd_truth = true(sim_params.data_reps, 1);
    H0_sig_grnd_truth = false(sim_params.data_reps, 1);

    % FALSE POSITIVE RATE
    % NaN value in SNR matrix means an invalid test and a false in H0 matrix
    % Which means H0 was wrongly rejected (false alarm)
    % This false value will contribute to the fa_cnt
    % This false value will also be included in the number of valid tests
    % Subtract number of NaNs from fa_cnt and total number of tests
    LLR_noise_NaNs = sum(isnan(LLR_noise),"all");
    
    % Calculate probability of false alarm
    % Wrongly rejected H0 when H0 was true
    % Ground truth is matrix of "trues" (accepting H0)
    % Therefore, interested when there is a "false" value
    FPR_vec = H0_noise_grnd_truth - H0_noise;
    FPR_cnt = sum(FPR_vec,"all") - LLR_noise_NaNs;
    valid_noise_tests = (sim_params.data_reps - LLR_noise_NaNs);
    FPR = FPR_cnt / valid_noise_tests;

    % FALSE NEGATIVE RATE
    % NaN value in SNR matrix means an invalid test and a true in H0 matrix
    % Which means H0 was wrongly accepted (missed detection)
    % This incorrect acceptance contributes to the md_cnt
    % This true value will also be included in the number of valid tests
    % Subtract number of NaNs from md_cnt and total number of tests
    LLR_sig_NaNs = sum(isnan(LLR_sig),"all");

    % Calculate probability of missed detection
    % Accepted H0 when H0 was not true
    % Ground truth is matrix of "falses"
    % Therefore, interested when there is a "true" value
    FNR_vec = H0_sig_grnd_truth - H0_sig;
    FNR_cnt = abs(sum(FNR_vec,"all"))- LLR_sig_NaNs;
    valid_sig_tests = (sim_params.data_reps - LLR_sig_NaNs);
    FNR = FNR_cnt / valid_sig_tests;

    % True positive rate is one minus false negative rate
    TPR = 1 - FNR;
    
end