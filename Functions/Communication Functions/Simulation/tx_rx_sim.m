function [rxSig_out, noise_sig_out, SNR_dB, txSig_delay_pwr, noise_pwr, BER_theo, rxSig_BER] = ...
    tx_rx_sim(x, z, sim_params)

    %% ALICE
    % Create system objects and generate signal to be transmitted

    % SIGNAL SOURCE
    % Outputs payload one sample or frame at a time
    % Create system object
    sig_src = dsp.SignalSource();
    % Set properties
    sig_src.Signal = sim_params.msg_bits; % payload
    sig_src.SamplesPerFrame = sim_params.payload_len;
    sig_src.SignalEndAction = 'Cyclic repetition';
    % Signal source
    txSig_bin = sig_src();
    
    % SCRAMBLE
    % Scramble input signal
    % Create system object
    scrambler = comm.Scrambler();
    % Set properties
    scrambler.CalculationBase = sim_params.scram_base;
    scrambler.Polynomial = sim_params.scram_poly;
    scrambler.InitialConditions = sim_params.scram_IC;
    % Scramble
    txSig_scram = scrambler(txSig_bin);

    % CONCATENATE
    % Append scrambled bit sequence to header
    txSig_concat = [sim_params.header; txSig_scram];
    
    % MODULATE
    % Modulation turns integers into complex numbers (I&Q values)
    % AD936x tracking algorithm expects a complex quadrature signal at the input
    % If statements to choose modulation scheme
    % Not necessary to create System Objects for QAM
    if sim_params.M == 2
        % BPSK MODULATOR
        % Create system object
        modulator = comm.BPSKModulator;
        % Set properties
        modulator.PhaseOffset = 0;
        modulator.OutputDataType = 'double';
        % Modulate using BPSK
        txSig_mod = modulator(txSig_concat);
        header_mod = modulator(sim_params.header);
    elseif sim_params.M == 4
        % QPSK MODULATOR
        % Create system object
        modulator = comm.QPSKModulator;
        % Set properties
        modulator.BitInput = true;
        modulator.PhaseOffset = pi/sim_params.M;
        modulator.SymbolMapping = 'Gray';
        modulator.OutputDataType = 'double';
        % Modulate using QPSK
        txSig_mod = modulator(txSig_concat);
        header_mod = modulator(sim_params.header);
    elseif sim_params.M == 8
        % Modulate using 8PSK
        txSig_mod = pskmod(txSig_concat, ...
                    sim_params.M, ...
                    pi/sim_params.M, ...
                    'gray', ...
                    'InputType', 'bit', ...
                    'OutputDataType', 'double');
        header_mod = pskmod(sim_params.header, ...
                    sim_params.M, ...
                    pi/sim_params.M, ...
                    'gray', ...
                    'InputType', 'bit', ...
                    'OutputDataType', 'double');
    elseif sim_params.M == 16 || sim_params.M == 32 || sim_params.M == 64
        % Modulate using QAM
        txSig_mod = qammod(txSig_concat, ...
            sim_params.M, ...
            'gray', ...
            'InputType','bit', ...
            'OutputDataType', 'double');
        header_mod = qammod(sim_params.header, ...
            sim_params.M, ...
            'gray', ...
            'InputType','bit', ...
            'OutputDataType', 'double');
    else
        disp('\nPlease input a valid modulation scheme.\n')
    end

    % SPREAD
    % Spread the modulated signal
    % Generate vector with "chip" repeated values for each modulated symbol
    txSig_rep = repelem(txSig_mod, sim_params.chip_no(z));
    header_rep = repelem(header_mod, sim_params.chip_no(z));
    % Generate chip sequence (comm.PNSequence?)
    chip_seq = 2*randi([0 1], sim_params.chip_no(z), 1) - 1;
    % Generate spreading sequences
    % Repeat the chip sequence for each symbol
    spr_seq = repmat(chip_seq, (sim_params.header_len + sim_params.payload_len)/sim_params.bps, 1);
    spr_seq_head = repmat(chip_seq, sim_params.header_len/sim_params.bps, 1);
    % Multiply extended vector of complex numbers (I&Q) with spreading sequence
    txSig_spread = txSig_rep.*spr_seq;
    header_spread = header_rep.*spr_seq_head;
    
    % RAISED COSINE FILTER
    % Pulse shaping by interpolating signal using raised-cosine FIR filter
    % Create system object
    tx_filter = comm.RaisedCosineTransmitFilter();
    % Set properties
    tx_filter.Shape = sim_params.pulse_shape;
    tx_filter.RolloffFactor = sim_params.alpha; % Excess bandwidth
    tx_filter.FilterSpanInSymbols = 10;
    tx_filter.OutputSamplesPerSymbol = sim_params.interpol; % Upsampling factor
    % Filter
    txSig_shaped = tx_filter(txSig_spread);

    %% CHANNEL
    % Create system objects
    % Adds WGN, delay, and PFO to inputted signal

    % Calculate frame size (should be length(txSig_spread))
    % MATLAB example multiplied this by interpol for use in channel equations?
    frame_size = (sim_params.header_len + sim_params.payload_len)*sim_params.chip_no(z)/sim_params.bps;

    % PHASE & FREQUENCY OFFSET
    % Apply phase and frequency offsets to input signal
    % Create system object and set properties
    pfo = comm.PhaseFrequencyOffset( ...
        'PhaseOffset',              sim_params.phase_offset, ...
        'FrequencyOffset',          sim_params.freq_offset, ...
        'SampleRate',               sim_params.f_s);

    % DELAY
    % Apply phase and frequency offsets to input signal
    % Create system object and set properties
    % Maximum delay must be less than 65535
    var_delay = dsp.VariableFractionalDelay( ...
        'MaximumDelay',             (frame_size*sim_params.interpol)/sim_params.chip_no(z), ...
        'OutputDataType',           'Same as first input');
    % Other relevant parameters
    delay_stp_size = 0.0125*(sim_params.interpol);
    delay_max = 2*(sim_params.interpol);
    delay_min = 0.1;
    
    %% BOB
    % Create system objects

    % INITIALISE
    
    % Preallocate despread signal vector
    rxSig_dspr = zeros((sim_params.header_len + sim_params.payload_len)/sim_params.bps, 1);
    % Parameters for coarse frequency offset estimate (PLL)
    mean_freq_offset = 0;
    p_cnt = 0;
    
    % Since we only calculate BER on message part, 000s are not...
    % ...necessary here, they are just place-holder
    BER_mask = zeros(sim_params.no_of_msgs*length(sim_params.msg)*7, 1);
    for i = 1 : sim_params.no_of_msgs
        BER_mask((i-1)*length(sim_params.msg)*7 + (1:length(sim_params.msg)*7)) = ...
        (i-1)*sim_params.msg_len*7 + (1:length(sim_params.msg)*7);
    end

    % AUTOMATIC GAIN CONTROL
    % Adaptively adjusts gain to achieve constant signal level at output
    % Useful for QAM or FM
    % Create system object
    agc = comm.AGC();
    % Set properties
    agc.DesiredOutputPower = 2;
    agc.AveragingLength = 50;
    agc.MaxPowerGain = 20;
    
    % RAISED COSINE FILTER
    % Pulse shaping by interpolating signal using raised-cosine FIR filter
    % Create system object
    rx_filter = comm.RaisedCosineReceiveFilter();
    % Set properties
    rx_filter.RolloffFactor = sim_params.alpha; % Excess bandwidth
    rx_filter.FilterSpanInSymbols = sim_params.filter_span;
    rx_filter.InputSamplesPerSymbol = sim_params.interpol;
    rx_filter.DecimationFactor = sim_params.decim_fact;
    
    % COARSE FREQUENCY ESTIMATION
    % Compensate for frequency offset of PAM, PSK, or QAM signal
    % Create system object
    crs_freq_est = comm.CoarseFrequencyCompensator();
    % Set properties
    crs_freq_est.Modulation = sim_params.mod_type;
    crs_freq_est.Algorithm = 'Correlation-based';
    crs_freq_est.MaximumFrequencyOffset = 6e3;
    crs_freq_est.SampleRate = sim_params.f_s/sim_params.decim_fact;
    
    % COARSE FREQUENCY COMPENSATION
    % Compensate for frequency offset of PAM, PSK, or QAM signal
    % Create system object
    crs_freq_comp = comm.PhaseFrequencyOffset();
    % Set properties
    crs_freq_comp.PhaseOffset = 0;
    crs_freq_comp.FrequencyOffsetSource = 'Input port';
    crs_freq_comp.SampleRate = sim_params.f_s/sim_params.decim_fact;
    
    % TIMING RECOVERY
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
    
    % FINE FREQUENCY COMPENSATION
    % Compensate for carrier frequency offset
    % Create system object
    fn_freq_comp = comm.CarrierSynchronizer();
    % Set properties
    fn_freq_comp.Modulation = sim_params.mod_type;
    fn_freq_comp.ModulationPhaseOffset = 'Auto';
    fn_freq_comp.SamplesPerSymbol = sim_params.interpol/sim_params.decim_fact; % Post-filter oversampling
    fn_freq_comp.DampingFactor = 1;
    fn_freq_comp.NormalizedLoopBandwidth = 0.01;
    
    % PREAMBLE DETECTION
    % Detect preamble in data
    % Create system object
    prb_det = comm.PreambleDetector();
    % Set properties
    % Feed preamble detector the modulated/spread preamble 'symbol' sequence
    prb_det.Input = 'Symbol';
    prb_det.Threshold = 20; % Play with this?
    prb_det.Preamble = header_spread;
    
    % FRAME SYNCHRONISATION
    % Converts fixed- or variable-size inputs...
    % ...into fixed-size output frames of specified length
    % Valid frame always starts with a preamble
    % Create system object
    frame_sync = comm.internal.examples.FrameSynchronizer();
    % Set properties
    frame_sync.OutputFrameLength = frame_size;
    frame_sync.PreambleLength = sim_params.header_len * sim_params.chip_no(z) / sim_params.bps;
    
    % DEMODULATION
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
    
    % DESCRAMBLE
    % Descramble input signal
    % Create system object
    descrambler = comm.Descrambler( ...
        sim_params.scram_base, ...
        sim_params.scram_poly, ...
        sim_params.scram_IC);
    
    % BIT ERROR RATE
    % Compute bit or symbol error rate of input data
    % Create system object
    bit_err_calc = comm.ErrorRate( ...
        'Samples', 'Custom', ...
        'CustomSamples', BER_mask);
    
    % CONSTELLATION DIAGRAM
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
    const_diagram.SymbolsToDisplay = sim_params.payload_len/sim_params.bps;
    % Display Configuration
    const_diagram.ColorFading = 0;
    const_diagram.XLimits = axis_limits;
    const_diagram.YLimits = axis_limits; 
    const_diagram.ShowLegend = true;
    const_diagram.ChannelNames = {'Channel 1'};
    % Reference Constellation
    const_diagram.ReferenceConstellation = ref_const;
    
    %% OPERATION    
    % Pass signal through channel continuously (enables PLL to work)

    % Number of frames transmitted depends on stop time
    frame_time = sim_params.T_sym * frame_size;
    num_frames = sim_params.stop_time_bob/frame_time;

    % Initialize BER vector in case no valid frame is found
    rxSig_BER = NaN(3,1);

    % Transmit multiple frames
    for count = 1 : num_frames
        %% CHANNEL
    
        % Signal undergoes phase/frequency offset
        txSig_rotate = pfo(txSig_shaped);
    
        % Calculate the delay 
        if strcmp(sim_params.delay_type,'Ramp')
            % Variable delay taking the form of a ramp
            delay = min(((count - 1) * delay_stp_size + delay_min), ...
                (frame_size - sim_params.filter_span) * sim_params.interpol);
        else
            % Variable delay taking the shape of a triangle
            idx = mod((count - 1), (2*delay_max/delay_stp_size));
            if idx <= delay_max/delay_stp_size
                delay = idx * delay_stp_size;
            else
                delay = (2 * delay_max) - (idx * delay_stp_size);
            end
        end
    
        % Delayed signal
        delay = cast(delay, 'like', real(txSig_shaped));
        txSig_delay = var_delay(txSig_rotate, delay);
    
        % Signal passing through AWGN channel
        % For a zero-mean signal, variance should be equal to power
        [txSig_corrupt, ~] = awgn(txSig_delay, sim_params.SNR_dB(x), 'measured', 'dB');

        %% RECEPTION
        
        % AGC of received signal
        %rxSig_agc = agc(txSig_corrupt);
                
        % SRRC Filtering
        rxSig_shaped = rx_filter(txSig_corrupt);
    
        % Average coarse frequency offset estimate
        % Enables carrier synchroniser to lock/converge
        [~,est_freq_offset] = crs_freq_est(rxSig_shaped);
        est_freq_offset = (est_freq_offset + p_cnt*mean_freq_offset)/(p_cnt + 1);
        % Update state
        p_cnt = p_cnt + 1;
        mean_freq_offset = est_freq_offset;
    
        % Coarse frequency compensation
        rxSig_crs_comp = crs_freq_comp(rxSig_shaped, -est_freq_offset);
        
        % Symbol timing recovery
        rxSig_sync = symbol_sync(rxSig_crs_comp);
        
        % Fine frequency compensation
        rxSig_fn_comp = fn_freq_comp(rxSig_sync);
    
        % Index of the peak corresponds to the end of the preamble
        [prb_idx, detmet] = prb_det(rxSig_fn_comp);
        
        % Frame synchronization
        [rxSig_frame_sync, frame_valid] = frame_sync(rxSig_fn_comp, prb_idx, detmet);
    
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
            
            %% PHASE AMBIGUITY COMPENSATION
            % Phase ambiguity estimation
            ph_ambiguity_est = round(angle(mean(conj(header_mod) .* ...
                rxSig_dspr(1:(sim_params.header_len/sim_params.bps))))*2/pi)/2*pi;    
            % Compensating for the phase ambiguity
            rxSig_res_ph = rxSig_dspr .* exp(-1i*ph_ambiguity_est);
    
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
       
            % Perform descrambling on payload portion only
            rxSig_descram = descrambler(rxSig_bin(sim_params.header_len + (1:sim_params.payload_len),1));
    
            % Recovering the message from the data
            char_set = int8(bit2int(rxSig_descram,7));
            if sim_params.print_bob == true
                fprintf('%s\n',char(char_set));
            end
    
            % Perform BER calculation only on message part
            rxSig_BER = bit_err_calc(sim_params.msg_bits, rxSig_descram);
        end
        
    end

    %% COMPUTE SNR

    % Obtain noise signal and it's power
    noise_sig = txSig_corrupt - txSig_delay;
    noise_pwr = mean(abs(noise_sig).^2);
    % Output signals of interest
    noise_sig_out = mat2cell(noise_sig, (frame_size*sim_params.interpol), 1); 
    rxSig_out = mat2cell(txSig_corrupt, (frame_size*sim_params.interpol), 1);
    % Calculate actual SNR
    txSig_delay_pwr = mean(abs(txSig_delay).^2);
    SNR_lin = txSig_delay_pwr/noise_pwr;
    SNR_dB = 10*log10(SNR_lin);
    
    %% THEORETICAL BER
    
    % Compare BER in terms of EbN0 because Eb depends on modulation scheme
    % Calculates probability of bit error based on Eb/No
    % Calculate Eb/No for each SNR value
    EbNo = SNR_lin*sim_params.BW*sim_params.T_sym/sim_params.bps;
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
    
    % Print results of interest
    fprintf('Signal-to-Noise ratio (dB) is = %f.\n', SNR_dB);
    fprintf('Theoretical bit error rate is = %f.\n', BER_theo);
    fprintf('Bit error rate is = %f.\n', rxSig_BER(1));
    fprintf('Number of detected errors = %d.\n', rxSig_BER(2));
    fprintf('Total number of compared samples = %d.\n', rxSig_BER(3));

    % Release system objects now that we are finished
    % Transmit objects
    release(sig_src)
    release(scrambler);
    release(tx_filter);
    if sim_params.M == 2 || sim_params.M == 4
        release(modulator);
    end
    % Channel objects
    release(pfo);
    release(var_delay);
    % Receive objects
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