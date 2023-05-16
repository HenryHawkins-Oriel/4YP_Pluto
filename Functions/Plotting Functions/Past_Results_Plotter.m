
% Load in whichever set of results you need
load('bob_results_.mat')

%% PLOT BER vs SNR
% Characterise Alice and Bob's channel
[SNR_dB_avrg, BER_avrg, BER_mat] = BER_vs_SNR(rxSig_BER, SNR_dB, 'median','sim', sim_params);
[SNR_dB_avrg, BER_avrg, BER_mat] = BER_vs_SNR(rxSig_BER, SNR_dB, 'median','real', real_params);
[SNR_dB_avrg, BER_theo_avrg] = BER_theo_vs_SNR(SNR_dB, BER_theo, 'mean','sim', sim_params);
[SNR_dB_avrg, BER_theo_avrg] = BER_theo_vs_SNR(SNR_dB, BER_theo, 'median','real', real_params);

%% PLOT BER vs Gain
% Characterise Alice and Bob's channel
[BER_avrg] = BER_vs_Gain(rxSig_BER, 'median', real_params);
% BER_theo vs Gain isn't really necessary/helpful
%[BER_theo_avrg] = BER_theo_vs_Gain(BER_theo, 'median','real', real_params);

%% PLOT LLR vs SNR
% Two line plots will allow us to find an accurate threshold in the middle
[A,B,C] = LLR_vs_SNR(SNR_dB, LLR_noise, LLR_sig, 'real', real_params);
[A,B,C] = LLR_vs_SNR(SNR_dB, LLR_noise, LLR_sig, 'sim', sim_params);

%% PLOT MDs or FAs vs SNR
[prob_md_vec, H0_sig] = MD_vs_SNR(B, SNR_dB, LLR_sig, 'mean','sim', sim_params);
[prob_md_vec, H0_sig] = MD_vs_SNR(B, SNR_dB, LLR_sig, 'mean','real', real_params);
[prob_fa_vec, H0_noise] = FA_vs_SNR(B, SNR_dB, LLR_noise, 'mean','sim', sim_params);
[prob_fa_vec, H0_noise] = FA_vs_SNR(B, SNR_dB, LLR_noise, 'mean','real', real_params);

%% PLOT BER vs MDs or FAs
[prob_md_vec, H0_sig] = BER_vs_MD(B, SNR_dB, rxSig_BER, LLR_sig, 'mean','sim', sim_params);
[prob_md_vec, H0_sig] = BER_vs_MD(B, SNR_dB, rxSig_BER, LLR_sig, 'mean','real', real_params);
[prob_fa_vec, H0_noise] = BER_vs_FA(B, SNR_dB, rxSig_BER, LLR_noise, 'mean','sim', sim_params);
[prob_fa_vec, H0_noise] = BER_vs_FA(B, SNR_dB, rxSig_BER, LLR_noise, 'mean','real', real_params); % ROC?

%% RECEIVER OPERATING CHARACTERISTIC
[TPR, FPR] = ROC(A, B, C, LLR_noise, LLR_sig, SNR_dB, 'sim', sim_params);
[TPR, FPR] = ROC(A, B, C, LLR_noise, LLR_sig, SNR_dB, 'real', real_params);

%% NOISE CHARACTERISTICS
% Bivariate histogram of received data
histogram2(real(noise_sig{1,1,1}),imag(noise_sig{1,1,1}),25,'FaceColor','flat');
xlabel('In-phase Component')
ylabel('Quadrature Component')
zlabel('Number of Samples')
colorbar

% Normplot
normplot(real(noise_sig{1,1,1}))
xlabel('In-phase Component')
normplot(imag(noise_sig{1,1,1}))
xlabel('Quadrature Component')

%% HEADER CHARACTERISTICS
[c,lags] = xcorr(sim_params.bip_barker);
stem(lags,c)
xlabel('Delay')
ylabel('Correlation')


