function [err_idx, num_errs, err_locs, LLR_sig, LLR_noise, rxSig_BER, BER_theo, SNR_dB, txSig_pwr] = ...
    replace_errs(err_log, LLR_sig, LLR_noise, rxSig_BER, BER_theo, SNR_dB, txSig_pwr)
    % Recognise when a test was invalid

    % Error statistics
    err_idx = find(~(cellfun(@isempty, {err_log.err_msg})));
    num_errs = length(err_idx);
    err_locs = transpose([err_log(err_idx).err_y;
                          err_log(err_idx).err_x;
                          err_log(err_idx).err_z]);

    % Ensure all results are NaN to indicate an invalid test
    if num_errs > 0
        for i = 1 : num_errs
             temp = err_locs(i,:);
             y = temp(1);
             x = temp(2);
             z = temp(3);
             % Ensure values are NaN
             LLR_sig(y,x,z) = NaN;
             LLR_noise(y,x,z) = NaN;
             rxSig_BER{y,x,z} = NaN(3,1);
             BER_theo(y,x,z) = NaN;
             SNR_dB(y,x,z) = NaN;
             txSig_pwr(y,x,z) = NaN;
        end
    end
    
end
