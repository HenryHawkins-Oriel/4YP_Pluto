function [SNR_dB, BER_theo, rxSig_BER, rxSig_out, txSig_pwr, current_time] = ...
rx_bob(z, frame_size, noise_pwr, freq_corr, sim_params)
    %% INITIALISE
    
    % Import spread or modulated header over from Alice's file
    load('Alice.mat','chip_seq','header_mod','header_spread')

    % Preallocate all vectors where applicable
    rxSig_dspr = zeros((sim_params.header_len + sim_params.payload_len)/sim_params.bps, 1);
    % Parameters for coarse frequency offset estimate
    mean_freq_offset = 0;
    p_cnt = 0;
    
    % Generate the target bits
    coder.extrinsic('sprintf');
    msg_set = zeros(sim_params.no_of_msgs*sim_params.msg_len, 1); 
    for msg_cnt = 0 : (sim_params.no_of_msgs - 1)
        msg_set(msg_cnt * sim_params.msg_len + (1 : sim_params.msg_len)) = ...
            sprintf('%s %03d\n', sim_params.msg, mod(msg_cnt, sim_params.no_of_msgs));
    end
    target_bits = int2bit(msg_set,7);
    
    % Since we only calculate BER on message part, 000s are not...
    % ...necessary here, they are just place-holder
    BER_mask = zeros(sim_params.no_of_msgs*length(sim_params.msg)*7, 1);
    for i = 1 : sim_params.no_of_msgs
        BER_mask((i-1)*length(sim_params.msg)*7 + (1:length(sim_params.msg)*7)) = ...
        (i-1)*sim_params.msg_len*7 + (1:length(sim_params.msg)*7);
    end
    
    %% AUTOMATIC GAIN CONTROL
    % Adaptively adjusts gain to achieve constant signal level at output
    % Useful for QAM or FM
    % Create system object
    agc = comm.AGC();
    % Set properties
    agc.DesiredOutputPower = 2;
    agc.AveragingLength = 50;
    agc.MaxPowerGain = 60;
    
    %% RAISED COSINE FILTER
    % Pulse shaping by interpolating signal using raised-cosine FIR filter
    % Create system object
    rx_filter = comm.RaisedCosineReceiveFilter();
    % Set properties
    rx_filter.RolloffFactor = sim_params.alpha; % Excess bandwidth
    rx_filter.FilterSpanInSymbols = 10;
    rx_filter.InputSamplesPerSymbol = sim_params.interpol;
    rx_filter.DecimationFactor = sim_params.decim_fact;
    
    %% COARSE FREQUENCY ESTIMATION
    % Compensate for frequency offset of PAM, PSK, or QAM signal
    
    % Create system object
    crs_freq_est = comm.CoarseFrequencyCompensator();
    % Set properties
    crs_freq_est.Modulation = sim_params.mod_type;
    crs_freq_est.Algorithm = 'Correlation-based';
    crs_freq_est.MaximumFrequencyOffset = 6e3;
    crs_freq_est.SampleRate = sim_params.f_s/sim_params.decim_fact;
    
    %% COARSE FREQUENCY COMPENSATION
    % Compensate for frequency offset of PAM, PSK, or QAM signal
    
    % Create system object
    crs_freq_comp = comm.PhaseFrequencyOffset();
    % Set properties
    crs_freq_comp.PhaseOffset = 0;
    crs_freq_comp.FrequencyOffsetSource = 'Input port';
    crs_freq_comp.SampleRate = sim_params.f_s/sim_params.decim_fact;
    
    %% TIMING RECOVERY
    % Correct symbol timing clock skew
    
    % Phase detector gain (Kp) for Timing Recovery PLL
    % Refer to "Digital Communications - A Discrete-Time Approach" by Michael Rice
    A = 1/sqrt(2);
    K = 1;
    if sim_params.M == 2
        K_p = 2.7*2*K*A^2; % 2.7
    elseif sim_params.M == 4
        K_p = 2*(2.7*2*K*A^2); % 5.4
    elseif sim_params.M == 8
        K_p = 2*(2*(2.7*2*K*A^2)); % 10.8
    elseif sim_params.M == 16
        K_p = 5.4; % Change
    elseif sim_params.M == 32
        K_p = 5.4; % Change
    elseif sim_params.M == 64
        K_p = 5.4; % Change
    else
        disp('\nPlease input a valid modulation scheme.\n')
    end
    
    % Create system object
    symbol_sync = comm.SymbolSynchronizer();
    % Set properties
    % Gardner operates at 2 samples/symbol for BPSK and QPSK
    % Gardner allows us to achieve timing synchronisation before...
    % ...carrier phase synchronisation
    symbol_sync.TimingErrorDetector = "Gardner (non-data-aided)"; 
    symbol_sync.SamplesPerSymbol = sim_params.interpol/sim_params.decim_fact; % Post-filter oversampling
    symbol_sync.DampingFactor = 1;
    symbol_sync.NormalizedLoopBandwidth = 0.01;
    symbol_sync.DetectorGain = K_p;
    
    %% FINE FREQUENCY COMPENSATION
    % Compensate for carrier frequency offset
    
    % Create system object
    fn_freq_comp = comm.CarrierSynchronizer();
    % Set properties
    fn_freq_comp.Modulation = sim_params.mod_type;
    fn_freq_comp.ModulationPhaseOffset = 'Auto';
    fn_freq_comp.SamplesPerSymbol = sim_params.interpol/sim_params.decim_fact; % Post-filter oversampling
    fn_freq_comp.DampingFactor = 1;
    fn_freq_comp.NormalizedLoopBandwidth = 0.01;
    
    %% PREAMBLE DETECTION
    % Detect preamble in data
    
    % Create system object
    prb_det = comm.PreambleDetector();
    % Set properties
    % Feed preamble detector the modulated/spread preamble 'symbol' sequence
    prb_det.Input = 'Symbol';
    prb_det.Threshold = 0.8; % Play with this?
    prb_det.Preamble = header_spread;
    
    %% FRAME SYNCHRONISATION
    % Converts fixed- or variable-size inputs...
    % ...into fixed-size output frames of specified length
    % Valid frame always starts with a preamble
    
    % Create system object
    frame_sync = comm.internal.examples.FrameSynchronizer();
    % Set properties
    frame_sync.OutputFrameLength = frame_size;
    frame_sync.PreambleLength = sim_params.header_len * sim_params.chip_no(z) / sim_params.bps;
    
    %% DEMODULATION
    
    % Create demodulator system objects
    % If statements choose modulation scheme
    if sim_params.M == 2
        % Create system object
        demodulator = comm.BPSKDemodulator;
        % Set properties
        demodulator.PhaseOffset = 0;
        demodulator.OutputDataType = 'double';
    elseif sim_params.M == 4
        % Create system object
        demodulator = comm.QPSKDemodulator;
        % Set properties
        demodulator.PhaseOffset = pi/sim_params.M;
        demodulator.BitOutput = true;
    end
    
    
    %% DESCRAMBLE
    % Descramble input signal
    % Create system object
    descrambler = comm.Descrambler( ...
        sim_params.scram_base, ...
        sim_params.scram_poly, ...
        sim_params.scram_IC);
    
    %% BIT ERROR RATE
    % Compute bit or symbol error rate of input data
    % Create system object
    bit_err_calc = comm.ErrorRate( ...
        'Samples', 'Custom', ...
        'CustomSamples', BER_mask);
    
    %% CONSTELLATION DIAGRAM
    
    % Constellation diagram parameters
    if sim_params.M == 2
        ref_const = pskmod(0:sim_params.M - 1, sim_params.M, 0);
        axis_limits = [-2 2];
    elseif sim_params.M == 4 || sim_params.M == 8
        ref_const = pskmod(0:sim_params.M - 1, sim_params.M, pi/sim_params.M);
        axis_limits = [-2 2];
    elseif sim_params.M == 16
        ref_const = qammod(0:sim_params.M - 1, sim_params.M);
        axis_limits = [-5 5];
    elseif sim_params.M == 32
        ref_const = qammod(0:sim_params.M - 1, sim_params.M);
        axis_limits = [-7 7];
    elseif sim_params.M == 64
        ref_const = qammod(0:sim_params.M - 1, sim_params.M);
        axis_limits = [-9 9];
    else
        disp('\nPlease input a valid modulation scheme.\n')
    end
    
    % Display and analyze input signals in IQ-plane
    % Create system object
    const_diagram = comm.ConstellationDiagram();
    % Set properties
    const_diagram.ShowReferenceConstellation = 1;
    const_diagram.ShowTrajectory = 0;
    const_diagram.EnableMeasurements = 0;
    const_diagram.NumInputPorts = 1;
    % Symbol configuration
    const_diagram.SampleOffset = 0;
    const_diagram.SamplesPerSymbol = sim_params.interpol/sim_params.decim_fact;
    const_diagram.SymbolsToDisplaySource = 'Property';
    const_diagram.SymbolsToDisplay = frame_size;
    % Display Configuration
    const_diagram.ColorFading = 0;
    const_diagram.XLimits = axis_limits;
    const_diagram.YLimits = axis_limits; 
    const_diagram.ShowLegend = false;
    const_diagram.ChannelNames = {'Channel 1'};
    % Reference Constellation
    const_diagram.ReferenceConstellation = ref_const;
    
    %% RECEPTION
    % Create Pluto receiver system object 
    rxPluto = sdrrx('Pluto');
    % Set radio properties
    rxPluto.RadioID = sim_params.Bob_Rx_ID;
    rxPluto.ShowAdvancedProperties = true;
    rxPluto.FrequencyCorrection = freq_corr; % ppm
    rxPluto.CenterFrequency = sim_params.f_c;
    rxPluto.BasebandSampleRate = sim_params.f_s;
    rxPluto.SamplesPerFrame = frame_size*sim_params.interpol;
    rxPluto.GainSource = 'Manual';
    rxPluto.OutputDataType = 'double';
    rxPluto.Gain = sim_params.rx_gain; % dB

    % Initialize variables
    current_time = 0;
    frame_time = sim_params.T_sym * frame_size;
    rxSig_BER = NaN(3,1);
    rxSig = complex(NaN((rxPluto.SamplesPerFrame),1));
    
    while current_time <  sim_params.stop_time_bob
        
        % Receive signal from the radio
        rxSig = rxPluto();

        % Update constellation diagram title
        const_diagram.Title = ("Raw " + sim_params.mod_type + " Signal at " + num2str(sim_params.f_c/1e6) + ' MHz');

        % AGC
        rxSig_agc = agc(rxSig);

        % Update constellation diagram title
        const_diagram.Title = (sim_params.mod_type + " Signal Post-AGC at " + ...
            num2str(sim_params.f_c/1e6) + ' MHz');
    
        % SRRC Filtering
        rxSig_shaped = rx_filter(rxSig_agc);

        % Update constellation diagram title
        const_diagram.Title = (sim_params.mod_type + " Signal Post-SRRC at " + ...
            num2str(sim_params.f_c/1e6) + ' MHz');
    
        % Average coarse frequency offset estimate
        % Enables carrier synchroniser to lock/converge
        [~,est_freq_offset] = crs_freq_est(rxSig_shaped);
        est_freq_offset = (est_freq_offset + p_cnt*mean_freq_offset)/(p_cnt + 1);
        % Update state
        p_cnt = p_cnt + 1;
        mean_freq_offset = est_freq_offset;
    
        % Coarse frequency compensation
        rxSig_crs_comp = crs_freq_comp(rxSig_shaped, -est_freq_offset);

        % Update constellation diagram title
        const_diagram.Title = (sim_params.mod_type + " Signal Post-Coarse-Frequency-Compensation at " + ...
            num2str(sim_params.f_c/1e6) + ' MHz');
        
        % Symbol timing recovery
        rxSig_sync = symbol_sync(rxSig_crs_comp);

        % Update constellation diagram title
        const_diagram.Title = (sim_params.mod_type + " Signal Post-Symbol-Timing-Recovery at " + ...
            num2str(sim_params.f_c/1e6) + ' MHz');
        
        % Fine frequency compensation
        rxSig_fn_comp = fn_freq_comp(rxSig_sync);

        % Update constellation diagram title
        const_diagram.Title = (sim_params.mod_type + " Signal Post-Fine-Frequency-Compensation at " + ...
            num2str(sim_params.f_c/1e6) + ' MHz');
    
        % Index of the peak corresponds to the end of the preamble
        [prb_idx, detmet] = prb_det(rxSig_fn_comp);
        
        % Frame synchronization
        [rxSig_frame_sync, frame_valid] = frame_sync(rxSig_fn_comp, prb_idx, detmet);

        % Update constellation diagram title
        const_diagram.Title = (sim_params.mod_type + " Signal Post-Frame-Synchronisation at " + ...
            num2str(sim_params.f_c/1e6) + ' MHz');
    
    
        if frame_valid   
            %% DESPREAD & ADD
            % For each symbol (with "chip" chips), multiply by chip sequence
            % Then, for each symbol, add all those chips together
            % Then divide the result by the number of chips
            % Can be done efficiently via vector multiplication
            % Adding averages out noise
            % Tranpose received signal for vector multiplication
            rxSig_frame_sync_T = transpose(rxSig_frame_sync);
            % Iterate over each symbol
                for j = 1 : (sim_params.header_len + sim_params.payload_len)/sim_params.bps
                    rxSig_dspr(j,1) = (rxSig_frame_sync_T(1,(((j-1)*sim_params.chip_no(z))+1:(j*sim_params.chip_no(z)))) * chip_seq)/sim_params.chip_no(z);
                end

            % Update constellation diagram title
            const_diagram.Title = (sim_params.mod_type + " Signal Post-De-Spreading at " + ...
                num2str(sim_params.f_c/1e6) + ' MHz');
    
            
            %% PHASE AMBIGUITY COMPENSATION
            % Phase ambiguity estimation
            ph_ambiguity_est = round(angle(mean(conj(header_mod) .* ...
                rxSig_dspr(1:(sim_params.header_len/sim_params.bps))))*2/pi)/2*pi;    
            % Compensating for the phase ambiguity
            rxSig_res_ph = rxSig_dspr .* exp(-1i*ph_ambiguity_est);

            % Update constellation diagram title
            const_diagram.Title = (sim_params.mod_type + " Signal Post-Phase-Ambiguity-Compensation at " + ...
                num2str(sim_params.f_c/1e6) + ' MHz');
    
            %% DEMODULATION
            % Demodulation returns complex values to integers, then to bits
            if sim_params.M == 2
                % Demodulate received signal using BPSK
                rxSig_demod = demodulator(rxSig_res_ph);
                rxSig_bin = int2bit(rxSig_demod,sim_params.bps);
            elseif sim_params.M == 4
                % Demodulate received signal using QPSK
                rxSig_bin = demodulator(rxSig_res_ph);
            elseif sim_params.M == 8
                % Demodulate received signal using 8PSK
                rxSig_bin = pskdemod(rxSig_res_ph, ...
                    sim_params.M, ...
                    pi/sim_params.M, ...
                    'gray', ...
                    'OutputType', 'bit');
            elseif sim_params.M == 16 || sim_params.M == 32 || sim_params.M == 64
            % 16-QAM, 32-QAM, or 64-QAM
                % Demodulate received signals using QAM
                rxSig_bin = qamdemod(rxSig_res_ph, ...
                    sim_params.M, ...
                    'gray', ...
                    'OutputType','bit');
            else
                disp('Please input a valid modulation scheme.')
            end

            % Update constellation diagram title
            const_diagram.Title = (sim_params.mod_type + " Signal Post-Demodulation at " + ...
                num2str(sim_params.f_c/1e6) + ' MHz');
       
            % Perform descrambling on payload portion only
            rxSig_descram = descrambler(rxSig_bin(sim_params.header_len + (1:sim_params.payload_len),1));

            % Update constellation diagram title
            const_diagram.Title = (sim_params.mod_type + " Signal Post-Descrambling at " + ...
                num2str(sim_params.f_c/1e6) + ' MHz');
    
            % Recovering the message from the data
            char_set = int8(bit2int(rxSig_descram,7));
            if sim_params.print_bob == true
                fprintf('%s\n',char(char_set));
            end
    
            % Perform BER calculation only on message part
            rxSig_BER = bit_err_calc(target_bits,rxSig_descram);

        end
    
        % Update simulation time
        current_time = current_time + frame_time;
        
    end

    %% COMPUTE SNR 
    % Calculate received signal power
    rxSig_pwr = mean(abs(rxSig).^2);
    txSig_pwr = (rxSig_pwr - noise_pwr);
    % Sometimes txSig_pwr comes out negative
    if txSig_pwr < 0
        % Disregard using NaN
        txSig_pwr = NaN;
        SNR_dB = NaN;
    else
        % Calculate average SNR
        SNR = txSig_pwr/noise_pwr;
        SNR_dB = 10*log10(SNR);
        % Calculate Eb/No for each SNR value
        EbNo = SNR*sim_params.BW*sim_params.T_sym/sim_params.bps;
    end
    
    %% THEORETICAL BER
    
    % Compare BER in terms of EbN0 because Eb depends on modulation scheme!
    % Calculates probability of bit error based on Eb/No
    if ~isnan(txSig_pwr)
        if sim_params.M == 2 || sim_params.M == 4 || sim_params.M == 8
            % BPSK, QPSK, or 8PSK
            BER_theo = (2/sim_params.bps)*qfunc(sqrt(EbNo*2*sim_params.bps*(sin(pi/sim_params.M))^2));
        elseif sim_params.M == 16 || sim_params.M == 32 || sim_params.M == 64
            % 16-QAM, 32-QAM, or 64-QAM
            BER_theo = (4/sim_params.bps)*((sqrt(sim_params.M) - 1)/sqrt(sim_params.M)) * ...
                qfunc(sqrt(EbNo*3*sim_params.bps/(sim_params.M - 1))); 
        else
            disp('Please input a valid modulation scheme.')
        end   
    else
        % Disregard using NaN
        BER_theo = NaN;
    end
    
    % Results of interest at end of signal acquisition
    rxSig_out = mat2cell(rxSig, rxPluto.SamplesPerFrame, 1);
    % Print
    if ~isnan(txSig_pwr)
        fprintf('Signal-to-Noise ratio (dB) is = %f.\n', SNR_dB);
        fprintf('Theoretical bit error rate is = %f.\n', BER_theo);
    else
        disp('Unable to calculate SNR and theoretical BER due to negative Eb/No.')
    end
    
    fprintf('Bit error rate is = %f.\n', rxSig_BER(1));
    fprintf('Number of detected errors = %d.\n', rxSig_BER(2));
    fprintf('Total number of compared samples = %d.\n', rxSig_BER(3));

    % Release system objects now that we are finished
    release(rxPluto);
    release(agc);
    release(rx_filter);
    release(crs_freq_est);
    release(crs_freq_comp);
    release(symbol_sync);
    release(fn_freq_comp);
    release(prb_det);
    release(frame_sync);
    release(descrambler);
    release(bit_err_calc);
    if sim_params.M == 2 || sim_params.M == 4
        release(demodulator);
    end

end