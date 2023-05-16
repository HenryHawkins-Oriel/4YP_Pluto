function RealParams = real_init
    % Parameters needed by all files

    % Pluto has been configured to operate using AD9364 firmware using PuTTY
    % New tuning range: 70 MHz to 6 GHz
    
    %% GENERAL PARAMETERS
    
    RealParams.M = 2;                                            % Modulation order (alphabet size)
    RealParams.bps = log2(RealParams.M);                         % Number of bits per symbol 
    RealParams.R_sym = 1e5;                                      % Symbol rate (Hz) (2e5 default)
    RealParams.T_sym = 1/RealParams.R_sym;                       % Symbol time (seconds)
    RealParams.interpol = 2;                                     % Interpolation factor
    RealParams.f_s = RealParams.R_sym * RealParams.interpol;     % Sample rate (scalar from 65,105 to 61.44e6)
    RealParams.f_c = 2.4e9;                                      % Centre or Carrier Frequency (915e6 or 2.4e9 or 5.2e9)
    RealParams.decim_fact = 1;                                   % Decimation factor
    RealParams.rx_gain = 70;                                     % Receiver gain (dB) ranges from -4 to 71 (must be 60 for 5.2 GHz)
    RealParams.radio_dist = 45e-2;                               % Measured distance between radios
    RealParams.tx_pwr = -23;                                     % Transmit power (dB) (7dBm)
        
    %% SCRAMBLER PARAMETERS
    
    RealParams.scram_base = 2;
    RealParams.scram_poly = [1 1 1 0 1];
    RealParams.scram_IC = [0 0 0 0];
    
    %% RAISED COSINE FILTER PARAMETERS

    RealParams.pulse_shape = 'Square root';         % Pulse shape for Raised Cosine (RC) Filters
    RealParams.filter_span = 10;                    % RC Filter span (in symbols)
    RealParams.alpha = 0.5;                         % Roll-off factor for SRRCs
    
    % Bandwidth calculation true for MPSK and MQAM
    if strcmp(RealParams.pulse_shape,'Square root')
        B = 0.5*(1 + RealParams.alpha);
        RealParams.BW = 2*B*RealParams.R_sym;
    else
        disp('Square Root Raised Cosine Filter must be chosen.')
    end
    
    %% PREAMBLE

    % [BarkerCode*2 | 'Hello world 000\n' | 'Hello world 001\n' ...];
    
    % Create a preamble using bipolar Barker code
    % Length must be an even number
    RealParams.bip_barker = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1]'; % Truncate +1 from end
    % Convert bipolar to unipolar (binary)
    RealParams.uni_barker = ((RealParams.bip_barker + 1)/2);
    % Send Barker code on both IQ components of QPSK modulated symbols
    % Achieved by repeating Barker code bits twice before modulating
    RealParams.header = repelem(RealParams.uni_barker,2);
    RealParams.header_len = length(RealParams.header);
    
    %% HELLO WORLD MESSAGES
    % 'Hello world 00000\n'...
    RealParams.msg = 'Hello world';
    RealParams.msg_len = length(RealParams.msg) + 5;
    % Number of messages in a frame
    % Ensure frame size is more than 3660 (better performance)
    RealParams.no_of_msgs = 33 * RealParams.bps;
    % 7 bits per characters (ASCII)
    RealParams.payload_len = RealParams.no_of_msgs * RealParams.msg_len * 7;
    
    % Converting characters to bits
    msg_set = zeros(RealParams.no_of_msgs * RealParams.msg_len, 1); 
    for msg_cnt = 0 : (RealParams.no_of_msgs - 1)
        msg_set(msg_cnt*RealParams.msg_len + (1:RealParams.msg_len)) = ...
        sprintf('%s %03d\n', RealParams.msg, msg_cnt);
    end
    bits = de2bi(msg_set,7,'left-msb')';
    RealParams.msg_bits = bits(:);

    %% MODULATION PARAMETERS

    % If statements to choose modulation scheme
    if RealParams.M == 2
        RealParams.mod_type = 'BPSK';
    elseif RealParams.M == 4
        RealParams.mod_type = 'QPSK';
    elseif RealParams.M == 8
        RealParams.mod_type = '8PSK';
    elseif RealParams.M == 16 || RealParams.M == 32 || RealParams.M == 64
        RealParams.mod_type = 'QAM';
    else
        disp('Please input a valid modulation scheme.')
    end

    %% STOP TIMES
    RealParams.ab_sync_time = 6;
    RealParams.stop_time_noise = 8; % Runs for +1s
    RealParams.stop_time_bob = 15; % Runs for +1s
    RealParams.stop_time_willie = 2; % Runs for +1s
    RealParams.stop_time_alice = (RealParams.stop_time_bob + 1) + 8;

    %% PRINT
    % Printing is nice, but it slows down the program and isn't essential
    RealParams.print_noise = false;
    RealParams.print_bob = false;
    RealParams.print_willie = false;

    %% FAIL SAFE
    % Fail-safe for while loops that synchronise Alice and Bob
    RealParams.count_max = 1000;

    %% RADIO IDs

    % Specify Radio IDs
    % 010 (NAB) for v4 (sw)
    % 000 (ABW) for v5 (dw) - use whichever is available
    RealParams.Alice_ID = 'usb:1';
    RealParams.Bob_Noise_ID = 'usb:0';
    RealParams.Bob_Rx_ID = 'usb:0';
    %SimParams.Willie_Noise_ID = 'usb:0';
    %SimParams.Willie_Rx_ID = 'usb:0';

    %% TRANSMITTER GAIN & SPREADING LOOP
    % Transmitter gain ranges from -89.75 to 0 dB with a resolution of 0.25 dB
    % Adalm Pluto transmitter power at 0 dB? Frequency dependent
    RealParams.tx_gain = -80:5:0; % Test ideal range (LLR 5 values)
    RealParams.data_reps = 15; % The more, the better
    RealParams.chip_no = [1,16,64]; % (1:1:8).^2 % Max chips per bit = 74


