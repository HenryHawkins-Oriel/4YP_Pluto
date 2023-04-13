function [noise_sig_out, noise_pwr, SNR_dB_theo, noise_fig, noise_flr, current_time] = ...
    rx_noise_sim(x, frame_size, freq_corr, radio_ID, sim_params)
    % Provides an estimate of the noise along our channel

    %% RECEPTION
    
    % Create Pluto receiver system object 
    rxPluto = sdrrx('Pluto');
    % Set radio properties
    rxPluto.RadioID = radio_ID;
    rxPluto.ShowAdvancedProperties = true;
    rxPluto.FrequencyCorrection = freq_corr; % ppm
    rxPluto.CenterFrequency = sim_params.f_c;
    rxPluto.BasebandSampleRate = sim_params.f_s;
    rxPluto.SamplesPerFrame = frame_size*sim_params.interpol;
    rxPluto.GainSource = 'Manual';
    rxPluto.OutputDataType = 'double';
    rxPluto.Gain = sim_params.rx_gain; % dB

    % Initialise current time
    current_time = 0;
    frame_time = sim_params.T_sym * frame_size;
    
    while current_time <  sim_params.stop_time_noise
        
        % Receive signal from the radio
        noise_sig = rxPluto();

        % Print noise samples
        if sim_params.print_noise == true
            disp(noise_sig);
        end
    
        % Update simulation time
        current_time = current_time + frame_time;
        
    end

    % Release radio system object so Alice can use this Radio ID
    release(rxPluto);

    % Noise signal outputted as cell
    noise_sig_out = mat2cell(noise_sig, rxPluto.SamplesPerFrame, 1);
    
    % Calculate noise power
    % Identical to L2 norm
    noise_pwr = mean(abs(noise_sig).^2);

    % Theoretical SNR (dB)
    % Use path loss equation and calculated noise power
    SNR_dB_theo = sim_params.tx_pwr - 10*log10(noise_pwr) ...
        + sim_params.tx_gain(x) + sim_params.rx_gain ...
        + 20*log10((3e8)/(4*pi*sim_params.f_c)) ...
        - 20*log10(sim_params.radio_dist);

    % Noise Figure (NF)
    % Does this match up with datasheet?
    noise_fig = 10*log10(noise_pwr) - 10*log10(1.38e-23 * 290 * sim_params.BW) ...
        - (sim_params.tx_gain(x) + sim_params.rx_gain);

    % Minimum Discernable Signal (A.K.A. Noise Floor)
    % Do not confuse with Minimum *Detectable* Signal (MDS)
    noise_flr = 10*log10(1.38e-23 * 290 * sim_params.BW) + noise_fig;
    %MDS = noise_flr + SNR_dB_req;

    % LILLIEFORS TEST
    % Quantify normality of data
    % Null hypothesis is that data is normally distributed
    % Result is 1 if the test rejects the null hypothesis
    % Small values of p cast doubt on the validity of the null hypothesis
    %H0_lil = lillietest(real(noise_sig));

    % MAXIMUM LIKELIHOOD ESTIMATION (MLE)
    % MLE function estimates mean and standard deviation for received data
    %[mu_sig, conf_ints] = mle(real(noise_sig),'Distribution','Normal','Alpha',0.05);
    % Isolate mu and sigma estimates into variables
    %mu = mu_and_sig(1);
    %sig = mu_and_sig(2);
    
end
