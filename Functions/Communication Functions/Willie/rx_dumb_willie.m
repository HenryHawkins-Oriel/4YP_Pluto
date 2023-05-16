function [LLR, txSig_pwr, willie_sig_out, current_time] = ...
rx_dumb_willie(frame_size, noise_pwr, sim_params)
    % Listens for a random amount of time, at random time instances
    % Calculates an LLR value

    %% RECEPTION
    
    % Generate random number of samples to listen for
    n = frame_size * sim_params.interpol;
    % Frame sizes less than 3660 can yield poor performance
    % Frame size must be even in "single channel" configuration
    num_samples = 2*randi([1830 n], 1, 1);
    
    % Create Pluto receiver system object 
    rxPluto = sdrrx('Pluto');
    % Set radio properties
    rxPluto.RadioID = sim_params.Willie_Rx_ID;
    rxPluto.CenterFrequency = sim_params.f_c;
    rxPluto.BasebandSampleRate = sim_params.f_s;
    rxPluto.SamplesPerFrame = num_samples;
    rxPluto.GainSource = 'Manual';
    rxPluto.OutputDataType = 'double';
    rxPluto.Gain = sim_params.rx_gain; % dB

    % Initialise current time
    current_time = 0;
    frame_time = sim_params.T_sym * num_samples;
    
    while current_time <  sim_params.stop_time_willie
        
        % Receive signal from the radio
        willie_sig = rxPluto();

        % Print noise samples
        if sim_params.print_willie == true
            disp(willie_sig);
        end
    
        % Update simulation time
        current_time = current_time + frame_time;
        
    end

    % Convert received signal to cell for output
    willie_sig_out = mat2cell(willie_sig, rxPluto.SamplesPerFrame, 1);
    % Calculate received signal power
    willie_sig_pwr = real(mean(conj(willie_sig).*willie_sig,1));
    txSig_pwr = (willie_sig_pwr - noise_pwr);

    % Disregard using NaN
    if txSig_pwr > 0
        LLR = LLR_calc(willie_sig, noise_pwr, txSig_pwr, sim_params);
    else
        txSig_pwr = NaN; % Limitation if we listen to noise
        LLR = NaN;
    end

    % Release radio system object so this Radio ID can be used again
    release(rxPluto);
    
end
