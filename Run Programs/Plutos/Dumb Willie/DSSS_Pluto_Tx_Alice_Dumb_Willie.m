% Program that transmits desired message repeatedly to Bob

% Make sure you've opened MATLAB with administrator privileges

clear

% Store and check information regarding attached radios
connected_radios_tx = findPlutoRadio;
if length(connected_radios_tx) == 3
    %% INITIALISE
    % Transmitter parameter structure
    sim_params = sim_init;
    % Preallocate output vectors
    signal_template = cell(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no));
    frame_size = NaN(1, length(sim_params.chip_no));
    time_template = zeros(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no), 'uint64');
    t_start_tx = time_template;
    t_end_tx = time_template;
    
    %% SYNCHRONISE
    
    % Reset stopwatch variable
    % Helps if previous run was cut short by error
    current_time_alice = 0;
    save('alice_time.mat','current_time_alice')
    
    % Synchronise Alice, Bob, and Willie's programs
    % Could go without this but it makes things easier
    % Can't load files inside a while loop
    load('pause.mat','time_elap')
    for count = 1 : sim_params.count_max
       if time_elap > sim_params.ab_sync_time
         break;
       end
       pause(sim_params.ab_sync_time/11) % Be careful with how often you load
       load('pause.mat','time_elap')
       disp(time_elap)
    end
    
    % Wait for Willie to stop listening to noise
    % After this, Alice and Bob's programs only depend on each other
    disp('Pause to wait for Willie to finish listening to noise.')
    load('noise_time_willie.mat','current_time_noise_willie');
    for count = 1 : sim_params.count_max
       if current_time_noise_willie ~= 0
         disp("Willie's first and only noise reception has ended")
         break; % terminate execution of this smaller for-loop
       end
       pause(0.5) % Be careful with loading too frequently
       load('noise_time_willie.mat','current_time_noise_willie')
    end
    
    %% ITERATE
    % Iterate over different transmitter gain values
    
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
            txSig = signal_template;
    
            for x = 1 : length(sim_params.tx_gain)
        
                % Pause to allow Bob to save 'noise_time.mat' for use in next sub-sub-loop
                pause(0.75)
            
                % For syncing up transmission/reception at the end of each loop
                % Need to reset current_time_alice to zero for Bob
                current_time_alice = 0;
                save('alice_time.mat','current_time_alice')
                
                % Wait for Bob to stop listening to noise
                disp('Pause to wait for Bob to finish listening to noise.')
                load('noise_time_bob.mat','current_time_noise_bob');
                for count = 1 : sim_params.count_max
                   if current_time_noise_bob ~= 0
                     disp(['Noise reception for Transmission ', num2str(y), '.', num2str(z), '.', num2str(x), ' has ended'])
                     break; % terminate execution of this smaller for-loop
                   end
                   pause(0.5) % Be careful with loading too frequently
                   load('noise_time_bob.mat','current_time_noise_bob')
                end
                
                % Start timer for one loop
                t_start_tx(y,x,z) = tic;
            
                % Transmit
                disp(['Transmission ', num2str(y), '.', num2str(z), '.', num2str(x), ' has started'])
                [txSig(y,x,z), current_time_alice] = tx_alice(x, z, frame_size(1,z), sim_params);
            
                if current_time_alice ~= 0
                    disp(['Transmission ', num2str(y), '.', num2str(z), '.', num2str(x), ' has ended'])
                    % Save current time so Bob knows when Alice is finished
                    save('alice_time.mat','current_time_alice')
                end   
                
                % Record time taken for one loop
                t_end_tx(y,x,z) = toc(t_start_tx(y,x,z));
            
                % Pause so Bob has time to realise Alice's transmission has ended
                pause(1)
        
            end
    
        end
    
    end
    
    % Save results to be safe, use -v7.3 MAT-file version to store large data
    save('alice_results_.mat', ...
        'txSig','t_end_tx', ...
        'sim_params',...
        '-v7.3')
    disp('**WARNING**: DO NOT RUN PROGRAM AGAIN WIHTOUT CHANGING RESULTS FILE NAME!')

else
    % End program if less than two radios are connected
    disp('Less than three radios are connected')
    return
end