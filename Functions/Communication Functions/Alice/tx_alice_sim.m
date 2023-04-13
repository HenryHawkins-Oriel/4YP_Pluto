function [txSig_out, chip_seq, header_mod, header_spread] = ...
tx_alice_sim(x, z, frame_size, sim_params)

    % Set rng function to default state so example produces repeatable results 
    rng default;
    
    %% DATA
    
    % SIGNAL SOURCE
    % Outputs payload one sample or frame at a time
    % Create system object
    sig_src = dsp.SignalSource();
    % Set properties
    sig_src.Signal = sim_params.msg_bits; % payload
    sig_src.SamplesPerFrame = sim_params.payload_len;
    sig_src.SignalEndAction = 'Cyclic repetition';
    
    txSig_bin = sig_src();
    
    %% SCRAMBLE
    % Scramble input signal
    % Create system object
    scrambler = comm.Scrambler();
    % Set properties
    scrambler.CalculationBase = sim_params.scram_base;
    scrambler.Polynomial = sim_params.scram_poly;
    scrambler.InitialConditions = sim_params.scram_IC;
    
    % Scramble
    txSig_scram = scrambler(txSig_bin);
    
    %% CONCATENATE
    % Append scrambled bit sequence to header
    txSig_concat = [sim_params.header; txSig_scram];
    
    %% MODULATE
    
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
        % Modulate
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
    
    %% SPREAD
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
    
    %% RAISED COSINE FILTER (TRANSMIT)
    
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

    %% TRANSMISSION

    % Calculate power of output signal for use in SNR calculation using CSI
    txSig_pwr_alice = mean(abs(txSig_shaped).^2);
    % The above is identical to both of the below
    %norm(txSig_shaped,2)^2/numel(txSig_shaped)
    %real(mean(conj(txSig_shaped).*txSig_shaped, 1))

    % Transmission Process (no need to transmit frames repeatedly)
    % Output signal as a cell for storage
    txSig_out = mat2cell(txSig_shaped, frame_size*sim_params.interpol, 1);
    

    % How do I include radio properties? Do I need to?
    BasebandSampleRate = sim_params.f_s;
    SamplesPerFrame = frame_size*sim_params.interpol;
    CenterFrequency = sim_params.f_c;
    Gain = sim_params.tx_gain(x); % dB
   
    
    % Release system objects now that we are finished
    release(scrambler);
    release(tx_filter);
    if sim_params.M == 2 || sim_params.M == 4
        release(modulator);
    end

end
