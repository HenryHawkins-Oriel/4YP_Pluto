function SimParams = sim_init
    % Parameters needed by all files

    % Pluto has been configured to operate using AD9364 firmware using PuTTY
    % New tuning range: 70 MHz to 6 GHz
    
    %% GENERAL PARAMETERS
    
    SimParams.M = 8;                                            % Modulation order (alphabet size)
    SimParams.bps = log2(SimParams.M);                          % Number of bits per symbol 
    SimParams.R_sym = 1e5;                                      % Symbol rate (Hz) (2e5 default)
    SimParams.T_sym = 1/SimParams.R_sym;                        % Symbol time (seconds)
    SimParams.interpol = 2;                                     % Interpolation factor
    SimParams.f_s = SimParams.R_sym * SimParams.interpol;       % Sample rate (scalar from 65,105 to 61.44e6)
    SimParams.decim_fact = 1;                                   % Decimation factor

    %% SCRAMBLER PARAMETERS
     
    SimParams.scram_base = 2;
    SimParams.scram_poly = [1 1 1 0 1];
    SimParams.scram_IC = [0 0 0 0];
    
    %% RAISED COSINE FILTER PARAMETERS

    SimParams.pulse_shape = 'Square root';         % Pulse shape for Raised Cosine (RC) Filters
    SimParams.filter_span = 10;                    % RC Filter span (in symbols)
    SimParams.alpha = 0.5;                         % Roll-off factor for SRRCs
    
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
    % Makes frame size roughly independent of bps
    % Ensure frame size is more than 3660 (better performance)
    SimParams.no_of_msgs = 33 * SimParams.bps;
    % 7 bits per characters (ASCII)
    SimParams.payload_len = SimParams.no_of_msgs * SimParams.msg_len * 7;
    
    % Converting characters to bits
    msg_set = zeros(SimParams.no_of_msgs * SimParams.msg_len, 1); 
    for msg_cnt = 0 : (SimParams.no_of_msgs - 1)
        msg_set(msg_cnt*SimParams.msg_len + (1:SimParams.msg_len)) = ...
            sprintf('%s %03d\n', SimParams.msg, mod(msg_cnt, SimParams.no_of_msgs)); % , msg_cnt)
    end
    SimParams.msg_bits = int2bit(msg_set, 7);

    %% CHANNEL PARAMETERS
    SimParams.phase_offset = 47;                % Degrees
    SimParams.freq_offset = 5e3;               % Hertz (Hz)
    SimParams.delay_type = 'Triangle';          % Type of delay for channel distortion

    %% MODULATION PARAMETERS

    % If statements to choose modulation scheme
    if SimParams.M == 2
        SimParams.mod_type = 'BPSK';
    elseif SimParams.M == 4
        SimParams.mod_type = 'QPSK';
    elseif SimParams.M == 8
        SimParams.mod_type = '8PSK';
    elseif SimParams.M == 16 || SimParams.M == 32 || SimParams.M == 64
        SimParams.mod_type = 'QAM';
    else
        disp('Please input a valid modulation scheme.')
    end

    %% STOP TIMES
    % This will dictate how many frames can be transmitted
    % Therefore what our maximum frame time (thus SF) can be, ...
    % ... as we need to fit in at least one frame
    % Beyond SF 134, BER drops to zero sharply (uninteresting)
    SimParams.stop_time_bob = 15; % Max SF of 403 (see paper working)

    %% PRINT
    % Printing is nice, but it slows down the program and isn't essential
    SimParams.print_bob = false;

    %% SNR & SPREADING LOOPB
    % Each iteration of simulation takes about 3.7s
    SimParams.SNR_dB = 100; % Test ideal range
    SimParams.data_reps = 2; % The more, the better
    SimParams.chip_no = 4; % (1:1:8).^2 (9 seems to perform badly)



