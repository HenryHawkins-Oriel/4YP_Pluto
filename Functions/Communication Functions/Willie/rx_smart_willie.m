function [null_hypo, LLR] = rx_smart_willie(sig, noise_pwr, txSig_pwr, sim_params)
    % Rejects or accepts the null hypothesis (no transmission, just noise)

    % Estimate threshold (experimentally)
    %threshold = log((qfuncinv(alice_params.prob_fa/2))^2);
    % Check suboptimal test

    if ~isnan(txSig_pwr)
        
        % Relevant test statistic parameters
        N = length(sig);
        noise_psd = noise_pwr/sim_params.BW; % W/Hz
        gamma_c = txSig_pwr*sim_params.T_sym;
    
        % Compute test statistic
        % LLR increases with SNR
        % Implementation of log(cosh(x)) that doesn't...
        % ...suffer from numberical errors
        x = (2*sqrt(txSig_pwr)/noise_psd) * abs(sig);
        lncosh_vec = log(0.5*ones(N,1)) + x + log(1 + exp(-2*x));
        lncosh_sum = sum(lncosh_vec);
        LLR = -(N*gamma_c) + lncosh_sum;

        % Compute threshold
        % Threshold is SNR dependent
        SNR = txSig_pwr/noise_pwr;
        SNR_dB = 10*log10(SNR);
        % Threshold equation found experimentally
        thresh = 10^(0.057*SNR_dB + 8.93);
    
        % Decision
        null_hypo = willie_decision(LLR_sig, thresh);
    else
        null_hypo = false;
        LLR = NaN;
    end
    
end
