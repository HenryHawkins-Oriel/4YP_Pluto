function SimParams = sim_init
    % Parameters needed by all files

    % Pluto has been configured to operate using AD9364 firmware using PuTTY
    % New tuning range: 70 MHz to 6 GHz
    
    %% GENERAL PARAMETERS
    
    SimParams.M = 2;                                            % Modulation order (alphabet size)
    SimParams.bps = log2(SimParams.M);                          % Number of bits per symbol 
    SimParams.R_sym = 1e5;                                      % Symbol rate (Hz) (2e5 default)
    SimParams.T_sym = 1/SimParams.R_sym;                        % Symbol time (seconds)
    SimParams.interpol = 2;                                     % Interpolation factor
    SimParams.f_s = SimParams.R_sym * SimParams.interpol;       % Sample rate (scalar from 65,105 to 61.44e6)
    SimParams.f_c = 5.2e9;                                      % Centre or Carrier Frequency (915e6 or 2.4e9 or 5.2e9)
    SimParams.decim_fact = 1;                                   % Decimation factor
    SimParams.pulse_shape = 'Square root';                      % Pulse shape for Raised Cosine Filter
    SimParams.alpha = 0.5;                                      % Roll-off factor for SRRCs
    SimParams.rx_gain = 60;                                     % Receiver gain (dB) ranges from -4 to 71 (must be 60 for 5.2 GHz)
    SimParams.radio_dist = 12e-2;                               % Measured distance between radios
    SimParams.tx_pwr = -23;                                     % Transmit power (dB) (7dBm)
        
    %% SCRAMBLER PARAMETERS
    
    SimParams.scram_base = 2;
    SimParams.scram_poly = [1 1 1 0 1];
    SimParams.scram_IC = [0 0 0 0];
    
    %% CALCULATE BANDWIDTH
    % Bandwidth calculation true for MPSK and MQAM
    if strcmp(SimParams.pulse_shape,'Square root')
        B = 0.5*(1 + SimParams.alpha);
        SimParams.BW = 2*B*SimParams.R_sym;
    else
        disp('Square Root Raised Cosine Filter must be chosen.')
    end
    
    %% PREAMBLE

    % [BarkerCode*2 | 'Hello world 000\n' | 'Hello world 001\n' ...];
    
    % Create a preamble using bipolar Barker code
    % Length must be an even number
    SimParams.bip_barker = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1]'; % Truncate +1 from end
    % Convert bipolar to unipolar (binary)
    SimParams.uni_barker = ((SimParams.bip_barker + 1)/2);
    % Send Barker code on both IQ components of QPSK modulated symbols
    % Achieved by repeating Barker code bits twice before modulating
    SimParams.header = repelem(SimParams.uni_barker,2);
    SimParams.header_len = length(SimParams.header);
    
    %% HELLO WORLD MESSAGES
    % 'Hello world 00000\n'...
    SimParams.msg = 'Hello world';
    SimParams.msg_len = length(SimParams.msg) + 5;
    % Number of messages in a frame
    % Ensure frame size is more than 3660 (better performance)
    if SimParams.M == 2
        SimParams.no_of_msgs = 33; % BPSK
    elseif SimParams.M == 4
        SimParams.no_of_msgs = 66; % QPSK
    elseif SimParams.M == 8
        SimParams.no_of_msgs = 99; % 8PSK
    else
        disp('Please input a valid modulation scheme.')
    end
    % 7 bits per characters (ASCII)
    SimParams.payload_len = SimParams.no_of_msgs * SimParams.msg_len * 7;
    
    % Converting characters to bits
    msg_set = zeros(SimParams.no_of_msgs * SimParams.msg_len, 1); 
    for msg_cnt = 0 : (SimParams.no_of_msgs - 1)
        msg_set(msg_cnt*SimParams.msg_len + (1:SimParams.msg_len)) = ...
        sprintf('%s %03d\n', SimParams.msg, msg_cnt);
    end
    bits = de2bi(msg_set,7,'left-msb')';
    SimParams.msg_bits = bits(:);

    %% MODULATION PARAMETERS

    % If statements to choose modulation scheme
    if SimParams.M == 2
        SimParams.mod_type = 'BPSK';
    elseif SimParams.M == 4
        SimParams.mod_type = 'QPSK';
    elseif SimParams.M == 8
        SimParams.mod_type = '8PSK';
    elseif SimParams.M == 16
        SimParams.mod_type = 'QAM';
    elseif SimParams.M == 32
        SimParams.mod_type = 'QAM';
    elseif SimParams.M == 64
        SimParams.mod_type = 'QAM';
    else
        disp('Please input a valid modulation scheme.')
    end

    %% STOP TIMES
    SimParams.ab_sync_time = 5;
    SimParams.stop_time_noise = 8; % Runs for +1s
    SimParams.stop_time_bob = 15; % Runs for +1s
    SimParams.stop_time_willie = 2; % Runs for +1s
    SimParams.stop_time_alice = (SimParams.stop_time_bob + 1) + 8;

    %% PRINT
    % Printing is nice, but it slows down the program and isn't essential
    SimParams.print_noise = false;
    SimParams.print_bob = false;
    SimParams.print_willie = false;

    %% FAIL SAFE
    % Fail-safe for while loops that synchronise Alice and Bob
    SimParams.count_max = 1000;

    %% RADIO IDs

    % Specify Radio IDs
    % 010 (NAB) for v4 (sw)
    % 000 (ABW) for v5 (dw) - use whichever is available
    SimParams.Alice_ID = 'usb:0';
    SimParams.Bob_Noise_ID = 'usb:0';
    SimParams.Bob_Rx_ID = 'usb:0';
    SimParams.Willie_Noise_ID = 'usb:0';
    SimParams.Willie_Rx_ID = 'usb:0';

    %% TRANSMITTER GAIN LOOP
    % Transmitter gain ranges from -89.75 to 0 dB with a resolution of 0.25 dB
    % Adalm Pluto transmitter power at 0 dB? Frequency dependent
    % Note that program struggles to get beyond 800 iterations
    % Unlikely to get past 500 iterations
    SimParams.tx_gain = -45:5:-20; % Test ideal range
    SimParams.data_reps = 1; % The more, the better
    SimParams.chip_no = 25; % (1:1:8).^2 % Max chips per bit = 74



