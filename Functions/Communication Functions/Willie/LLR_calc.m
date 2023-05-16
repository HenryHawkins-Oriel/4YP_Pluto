function LLR = LLR_calc(sig, noise_pwr, txSig_pwr, params)
    % Rejects or accepts the null hypothesis (no transmission, just noise)

    % Estimate threshold (experimentally)
    % Check suboptimal test?

    if ~isnan(txSig_pwr)
        % Relevant test statistic parameters
        N = length(sig);
        noise_psd = noise_pwr/params.BW; % W/Hz
        gamma_c = txSig_pwr*params.T_sym;
    
        % Compute test statistic
        % LLR increases with SNR
        % Implementation of log(cosh(x)) that doesn't...
        % ...suffer from numerical errors
        x = (2*sqrt(txSig_pwr)/noise_psd) * abs(sig);
        %lncosh_vec = x + log(0.5*(1 + exp(-2*x)));
        lncosh_sum = sum(x + log(0.5*(1 + exp(-2*x))));
        %lncosh_sum = sum(lncosh_vec);
        LLR = -(N*gamma_c) + lncosh_sum;
    else
        LLR = NaN;
    end
    
end
