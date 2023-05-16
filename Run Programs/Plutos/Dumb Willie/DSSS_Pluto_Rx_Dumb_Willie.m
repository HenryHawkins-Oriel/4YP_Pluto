% Program that runs a composite synchronous coherent detector

% Make sure you've opened MATLAB with administrator privileges
% Start Willie first

clear

% Store and check information regarding attached radios
connected_radios_rx = findPlutoRadio;
if length(connected_radios_rx) == 3
    %% INITIALISE
    % Transmitter parameter structure
    sim_params = sim_init;
    % Preallocate output vectors
    % Using NaN vectors helps if an experiment is cut short
    % Signal powers
    results_template = NaN(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no));
    txSig_pwr = results_template;
    % Willie's hypothesis tests
    willie_wait = results_template;
    LLR_noise = results_template;
    LLR_sig = results_template;
    % Preallocate signal vectors
    signal_template = cell(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no));
    frame_size = NaN(1, length(sim_params.chip_no));
    % Preallocate stopwatch vectors
    time_template = zeros(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no), 'uint64');
    t_start_willie = time_template;
    t_end_willie = time_template;
    
    %% SYNCHRONISE
    % Reset stopwatch variable
    time_elap = 0;
    save('pause.mat','time_elap')
    current_time_noise_willie = 0;
    save('noise_time_willie.mat','current_time_noise_willie')
    
    % Synchronise Alice, Bob, and Willie's programs
    count = 0;
    tic
    while (toc < sim_params.ab_sync_time) && (count < sim_params.count_max)
        pause(sim_params.ab_sync_time/10)
        time_elap = toc
        save('pause.mat','time_elap')
        count = count + 1;
    end
    
    %% NOISE RECEPTION
    % Willie listens to noise once at beginning of his program
    % Assumes noise channel is static
    
    disp("Obtaining first and only noise sample")
    [~, noise_pwr_willie, ~, ~, ~, current_time_noise_willie] = ...
        rx_noise(1, 5000, 0, sim_params.Willie_Noise_ID, sim_params);
    
    if current_time_noise_willie ~= 0
        disp("Willie's first and only noise reception has ended")
        % Save current time so Alice and Bob know when Willie is finished
        save('noise_time_willie.mat','current_time_noise_willie')
    end 
    
    % Pause so Alice and Bob have time to realise Willie has finished receiving noise
    pause(1.5)
    
    %% TEST
    % Iterate over different chip and receiver gain values, repeatedly
    % Listen at random time instances, for a random number of samples
    
    for y = 1 : sim_params.data_reps
    
        for z = 1 : length(sim_params.chip_no)
    
            % Spreading will dictate frame size (and thus frame time)
            % Can be an even positive integer from 2 to 16,777,216
            % Using values less than 3660 can yield poor performance
            frame_size(1,z) = (sim_params.header_len + sim_params.payload_len)*sim_params.chip_no(z)/sim_params.bps;
    
            % Preallocate 3D cell array of signal vectors
            % frame_size changes with each "page"
            signal_template(:,:,z) = repmat({complex(zeros(frame_size(1,z) * sim_params.interpol, 1), 0)}, ...
                sim_params.data_reps, length(sim_params.tx_gain));
            willie_sig = signal_template;
    
            % NB: No frequency correction
    
            for x = 1 : length(sim_params.tx_gain)
    
                % Wait for a random amount of time before listening
                % Whilst ensuring that Willie finishes listening before...
                % ...Alice is finished transmitting
                % Whole program lasts roughly 30s
                max_wait = 23 - sim_params.stop_time_willie;
                wait_gap_min = sim_params.stop_time_noise - 2;
                wait_gap_max = sim_params.stop_time_noise + 2; % Might catch beginning of A's signal
                % Willie cannot wait between 6-11 secs (transition period)
                willie_wait_temp = randi([1, wait_gap_max], 1, 1);
                wait_add = randi([(wait_gap_max - wait_gap_min), (max_wait - willie_wait_temp)], 1, 1);
                willie_wait_temp(willie_wait_temp > wait_gap_min) = willie_wait_temp(willie_wait_temp > wait_gap_min) + wait_add;
                willie_wait(y,x,z) = willie_wait_temp;
                disp(['Willie will listen in ', num2str(willie_wait(y,x,z)), ' seconds'])
                pause(willie_wait(y,x,z))
    
                load('alice_time.mat','current_time_alice');
                load('noise_time_bob.mat','current_time_noise_bob');
                
                if current_time_noise_bob ~= 0 && current_time_alice == 0
            
                    % Start timer for one noise loop
                    t_start_willie(y,x,z) = tic;
            
                    % Signal is being transmitted
                    disp(['Willie is listening to Transmission ', num2str(y), '.', num2str(z), '.', num2str(x)])
                    [LLR_sig(y,x,z), txSig_pwr(y,x,z), willie_sig(y,x,z)] = ...
                        rx_dumb_willie(frame_size(1,z), noise_pwr_willie, sim_params);
    
            
                    % Record time taken for one loop
                    t_end_willie(y,x,z) = toc(t_start_willie(y,x,z));
            
                else    
    
                    % Start timer for one noise loop
                    t_start_willie(y,x,z) = tic;
    
                    % Just listening to noise
                    disp(['Willie is listening to the noise before Transmission ', num2str(y), '.', num2str(z), '.', num2str(x)])
                    [LLR_noise(y,x,z), txSig_pwr(y,x,z), willie_sig(y,x,z)] = ...
                    rx_dumb_willie(frame_size(1,z), noise_pwr_willie, sim_params);
    
                    % Record time taken for one loop
                    t_end_willie(y,x,z) = toc(t_start_willie(y,x,z));
    
                end
    
                % Wait for Alice to stop transmitting
                % Equivalent to waiting for end of sub-sub-loop
                disp('Pause to wait for Alice to finish transmitting.')
                load('alice_time.mat','current_time_alice');
                for count = 1 : sim_params.count_max
                   if current_time_alice ~= 0
                     disp(['Transmission ', num2str(y), '.', num2str(z), '.', num2str(x),' has ended'])
                     break; % terminate execution of this smaller for-loop
                   end
                   pause(0.5) % Be careful with loading too frequently
                   load('alice_time.mat','current_time_alice')
                end
    
            end
    
        end
    
    end
    
    % Save results to be safe, use -v7.3 MAT-file version to store large data
    save('willie_results_.mat', ...
        'willie_sig','txSig_pwr','noise_pwr_willie', ...
        'LLR_noise','LLR_sig',...
        't_end_willie', ...
        'frame_size', ...
        'sim_params',...
        '-v7.3')
    disp('**WARNING**: DO NOT RUN PROGRAM AGAIN WIHTOUT CHANGING RESULTS FILE NAME!')
    
    % Verify probability of false alarm over multiple iterations
    % Produce receiver operating characteristic (ROC)
    % Obtain curves for different SNR values (tx_gain)
    % Obtain data points for different thresholds (prob_fa)
    
else
    % End program if less than three radios are connected
    disp('Less than three radios are connected')
    return
end
