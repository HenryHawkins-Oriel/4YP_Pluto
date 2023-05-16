function [txSig_out, current_time] = tx_alice(x, z, frame_size, real_params)
    
    % Set rng function to default state so example produces repeatable results 
    rng default;
    
    %% DATA
    
    % SIGNAL SOURCE
    % Outputs payload one sample or frame at a time
    % Create system object
    sig_src = dsp.SignalSource();
    % Set properties
    sig_src.Signal = real_params.msg_bits; % payload
    sig_src.SamplesPerFrame = real_params.payload_len;
    sig_src.SignalEndAction = 'Cyclic repetition';
    
    txSig_bin = sig_src();
    
    %% SCRAMBLE
    % Scramble input signal
    % Create system object
    scrambler = comm.Scrambler();
    % Set properties
    scrambler.CalculationBase = real_params.scram_base;
    scrambler.Polynomial = real_params.scram_poly;
    scrambler.InitialConditions = real_params.scram_IC;
    
    % Scramble
    txSig_scram = scrambler(txSig_bin);
    
    %% CONCATENATE
    % Append scrambled bit sequence to header
    txSig_concat = [real_params.header; txSig_scram];
    
    %% MODULATE
    
    % Modulation turns integers into complex numbers (I&Q values)
    % AD936x tracking algorithm expects a complex quadrature signal at the input
    
    % If statements to choose modulation scheme
    % Not necessary to create System Objects for QAM
    if real_params.M == 2
        % BPSK MODULATOR
        % Create system object
        modulator = comm.BPSKModulator;
        % Set properties
        modulator.PhaseOffset = 0;
        modulator.OutputDataType = 'double';
        % Modulate
        txSig_mod = modulator(txSig_concat);
        header_mod = modulator(real_params.header);
    elseif real_params.M == 4
        % QPSK MODULATOR
        % Create system object
        modulator = comm.QPSKModulator;
        % Set properties
        modulator.BitInput = true;
        modulator.PhaseOffset = pi/real_params.M;
        modulator.SymbolMapping = 'Gray';
        modulator.OutputDataType = 'double';
        % Modulate using QPSK
        txSig_mod = modulator(txSig_concat);
        header_mod = modulator(real_params.header);
    elseif real_params.M == 8
        % Modulate using 8PSK
        txSig_mod = pskmod(txSig_concat, ...
                    real_params.M, ...
                    pi/real_params.M, ...
                    'gray', ...
                    'InputType', 'bit', ...
                    'OutputDataType', 'double');
        header_mod = pskmod(real_params.header, ...
                    real_params.M, ...
                    pi/real_params.M, ...
                    'gray', ...
                    'InputType', 'bit', ...
                    'OutputDataType', 'double');
    elseif real_params.M == 16 || real_params.M == 32 || real_params.M == 64
        % Modulate using QAM
        txSig_mod = qammod(txSig_concat, real_params.M,'gray','InputType','bit');
        header_mod = qammod(real_params.header, real_params.M,'gray','InputType','bit');
    else
        disp('\nPlease input a valid modulation scheme.\n')
    end
    
    %% SPREAD
    % Generate vector with "chip" repeated values for each modulated symbol
    txSig_rep = repelem(txSig_mod, real_params.chip_no(z));
    header_rep = repelem(header_mod, real_params.chip_no(z));
    % Generate chip sequence (comm.PNSequence?)
    chip_seq = 2*randi([0 1], real_params.chip_no(z), 1) - 1;
    % Generate spreading sequences
    % Repeat the chip sequence for each symbol
    spr_seq = repmat(chip_seq, (real_params.header_len + real_params.payload_len)/real_params.bps, 1);
    spr_seq_head = repmat(chip_seq, real_params.header_len/real_params.bps, 1);
    % Multiply extended vector of complex numbers (I&Q) with spreading sequence
    txSig_spread = txSig_rep.*spr_seq;
    header_spread = header_rep.*spr_seq_head;
    
    %% RAISED COSINE FILTER (TRANSMIT)
    
    % Pulse shaping by interpolating signal using raised-cosine FIR filter
    % Create system object
    tx_filter = comm.RaisedCosineTransmitFilter();
    % Set properties
    tx_filter.Shape = real_params.pulse_shape;
    tx_filter.RolloffFactor = real_params.alpha; % Excess bandwidth
    tx_filter.FilterSpanInSymbols = real_params.filter_span;
    tx_filter.OutputSamplesPerSymbol = real_params.interpol; % Upsampling factor
    
    % Filter
    txSig_shaped = tx_filter(txSig_spread);
    
    %% TRANSMISSION
    
    % Create Pluto transmitter system object 
    txPluto = sdrtx('Pluto');
    % Set radio properties
    txPluto.RadioID = real_params.Alice_ID;
    txPluto.BasebandSampleRate = real_params.f_s;
    txPluto.SamplesPerFrame = frame_size*real_params.interpol;
    txPluto.CenterFrequency = real_params.f_c;
    txPluto.Gain = real_params.tx_gain(x); % dB
    
    % Output signal as a cell for storage
    txSig_out = mat2cell(txSig_shaped, txPluto.SamplesPerFrame, 1);

    % Calculate power of output signal for use in SNR calculation using CSI
    %txSig_pwr_alice = mean(abs(txSig_shaped).^2);
    % The above is identical to both of the below
    %norm(txSig_shaped,2)^2/numel(txSig_shaped)
    %real(mean(conj(txSig_shaped).*txSig_shaped, 1))
    
    % Save necessary parameters for use in Bob's file
    save('Alice.mat', ...
        'chip_seq', ...
        'header_mod', ...
        'header_spread')
    
    % Transmit frames repeatedly
    current_time = 0;
    frame_time = real_params.T_sym * frame_size;
        
    % Transmission Process
    while current_time < real_params.stop_time_alice
    
        % Data transmission
        txPluto(txSig_shaped);
    
        % Update simulation time
        current_time = current_time + frame_time;

    end
    
    % Release system objects now that we are finished transmitting
    release(txPluto);
    release(scrambler);
    release(tx_filter);
    if real_params.M == 2 || real_params.M == 4
        release(modulator);
    end

end
