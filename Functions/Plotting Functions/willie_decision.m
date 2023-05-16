function H0 = willie_decision(LLR, thresh)
    % Rejects or accepts null hypothesis bases on LLR being exceeding or 
    % falling below a designated threshold
    if LLR > thresh
        % Transmission detected (positive)
        % Reject null hypothesis
        H0 = false;
    else
        % Noise deteted (negative)
        % Accept null hypothesis
        H0 = true;
    end
    
end