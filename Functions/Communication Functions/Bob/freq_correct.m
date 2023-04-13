function freq_corr_factor = freq_correct(frame_size, sim_params)
    % Frequency Correction for ADALM-PLUTO Radio
    
    %% PARAMETER INITIALISATION

    f_ref = 80e3;
    num_samples = frame_size * sim_params.interpol;
    s3 = exp(1j*2*pi*f_ref*[0:10000-1]'/sim_params.f_s);  % 80 kHz
    s = 0.6*s3/max(abs(s3)); % Scale signal to avoid clipping in the time domain
    
    %% TRANSMITTER SETUP
    % Default value of 0 for FrequencyCorrection
    % Corresponds to the factory-calibrated condition

    % Create Pluto transmitter system object 
    tx = sdrtx('Pluto');
    % Set radio properties
    tx.RadioID = 'usb:1';
    tx.BasebandSampleRate = sim_params.f_s;
    tx.SamplesPerFrame = num_samples;
    tx.CenterFrequency = sim_params.f_c;
    tx.Gain = 0; % This has an impact on result past a certain point (poor SNR)
    tx.ShowAdvancedProperties = true;
    
    % Send tone at 80 kHz
    transmitRepeat(tx, s);
    
    %% RECEIVER SETUP
    % Default value of 0 for FrequencyCorrection
    % Corresponds to the factory-calibrated condition

    % Create Pluto receiver system object 
    rx = sdrrx('Pluto');
    % Set radio properties
    rx.RadioID = 'usb:0';
    rx.CenterFrequency = sim_params.f_c;
    rx.BasebandSampleRate = sim_params.f_s;
    rx.SamplesPerFrame = num_samples;
    rx.GainSource = 'Manual';
    rx.OutputDataType = 'double';
    rx.Gain = sim_params.rx_gain; % dB
    rx.ShowAdvancedProperties = true;

    % Receive signal
    rxSig = rx();
    
    % Find the tone that corresponds to the 80 kHz transmitted tone
    y = fftshift(abs(fft(rxSig)));
    [~, idx] = findpeaks(y,'MinPeakProminence',max(0.5*y));
    f_received = (max(idx) - num_samples/2 - 1)/num_samples * sim_params.f_s; % Hz
    
    %% ESTIMATION
    % Estimate value of FrequencyCorrection
    freq_corr_factor = (f_received - f_ref) / (sim_params.f_c + f_ref) * 1e6;

    % Release radio system objects for later use
    release(tx);
    release(rx);

end 