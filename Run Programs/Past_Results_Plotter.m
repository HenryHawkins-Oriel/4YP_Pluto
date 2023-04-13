
load('bob_results_.mat')

% Eventually, save figures

%% PLOT BER vs SNR
% Characterise Alice and Bob's channel
[SNR_dB_avrg_sort, BER_avrg_sort] = BER_vs_SNR(rxSig_BER, SNR_dB, 'median', sim_params);
[~, BER_theo_avrg_sort] = BER_theo_vs_SNR(SNR_dB, BER_theo, 'median', sim_params);

%% PLOT BER vs Gain
% Characterise Alice and Bob's channel
[BER_avrg] = BER_vs_Gain(rxSig_BER, 'median', sim_params);
[BER_theo_avrg] = BER_theo_vs_Gain(BER_theo, 'median', sim_params);

%% PLOT LLR vs SNR
% Two line plots will allow us to find an accurate threshold in the middle
[A,B,C] = LLR_vs_SNR(SNR_dB, LLR_noise, LLR_sig, 'median', sim_params);

%% PLOT MDs or FAs vs SNR
[prob_md_vec, H0_sig] = MD_vs_SNR(B, SNR_dB, LLR_sig, 'mean', sim_params);
[prob_fa_vec, H0_noise] = FA_vs_SNR(B, SNR_dB, LLR_noise, 'mean', sim_params);

%% PLOT BER vs MDs or FAs
[prob_md_vec, H0_sig] = BER_vs_MD(B, SNR_dB, rxSig_BER, LLR_sig, 'mean', sim_params);
[prob_fa_vec, H0_noise] = BER_vs_FA(B, SNR_dB, rxSig_BER, LLR_noise, 'mean', sim_params);

%% RECEIVER OPERATING CHARACTERISTIC
[TPR, FPR, b] = ROC(A, B, C, LLR_noise, LLR_sig, SNR_dB, sim_params);

%% NOISE CHARACTERISTICS
% Bivariate histogram of received data
histogram2(real(noise_sig{1,1,1}),imag(noise_sig{1,1,1}),25,'FaceColor','flat');
colorbar

% noise_psd = pwelch(noise_sig)
%any()
%normplot(real(noise_sig{1,1,1}))
%normplot(imag(noise_sig{1,1,1}))
%[c,lags] = xcorr(noise_sig{1,1,1});
%stem(lags,c)


