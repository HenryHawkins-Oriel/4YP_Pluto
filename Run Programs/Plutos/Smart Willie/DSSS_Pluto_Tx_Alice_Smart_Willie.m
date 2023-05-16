% Program that transmits desired message repeatedly to Bob

% Make sure you've opened MATLAB with administrator privileges

clear

% Store and check information regarding attached radios
connected_radios_tx = findPlutoRadio;
if length(connected_radios_tx) == 2
    %% INITIALISE
    % Transmitter parameter structure
    real_params = real_init;
    % Preallocate output vectors
    signal_template = cell(real_params.data_reps, ...
                           length(real_params.tx_gain), ...
                           length(real_params.chip_no));
    results_template = NaN(real_params.data_reps, ...
                           length(real_params.tx_gain), ...
                           length(real_params.chip_no));
    frame_size = NaN(1, length(real_params.chip_no));
    % Alice's error log
    program_it = length(real_params.tx_gain) * ...
                        real_params.data_reps * ...
                        length(real_params.chip_no);
    alice_err_log = struct('err_msg', cell(1, program_it), ...
                           'err_func1', cell(1, program_it), ...
                           'err_func2', cell(1, program_it), ...
                           'err_y', cell(1, program_it), ...
                           'err_x', cell(1, program_it), ...
                           'err_z', cell(1, program_it));

    % Stopwatch vectors
    time_template = zeros(real_params.data_reps, length(real_params.tx_gain), length(real_params.chip_no), 'uint64');
    t_start_tx = time_template;
    t_end_tx = time_template;

    % Vectors and cell arrays that can only be defined with (varying) frame size
    for z = 1 : length(real_params.chip_no)

        % Spreading will dictate frame size (and thus frame time)
        % Can be an even positive integer from 2 to 16,777,216
        % Using values less than 3660 can yield poor performance
        frame_size(1,z) = (real_params.header_len + real_params.payload_len)*real_params.chip_no(z)/real_params.bps;
    
        % Preallocate 3D cell array of signal vectors
        % frame_size changes with each "page"
        signal_template(:,:,z) = repmat({complex(zeros(frame_size(1,z) * real_params.interpol, 1), 0)}, ...
            real_params.data_reps, length(real_params.tx_gain));
        txSig = signal_template;

    end
    
    %% SYNCHRONISE
    
    % Reset stopwatch variables
    % Helps if previous run was cut short by error
    current_time_tx_alice = 0;
    save('alice_time.mat','current_time_tx_alice')
    time_elap = 0;
    save('pause.mat','time_elap')
    
    % Synchronise Alice and Bob's programs
    count = 0;
    tic
    while (toc < real_params.ab_sync_time) && (count < real_params.count_max)
        pause(0.15) % Can save pretty often now (thanks to try/catch)
        time_elap = toc
        try
            save('pause.mat','time_elap')
        catch
            continue
        end
        count = count + 1;
    end
    
    %% ITERATE
    % Iterate over different transmitter gain values
    
    for y = 1 : real_params.data_reps
    
        for z = 1 : length(real_params.chip_no)
    
            for x = 1 : length(real_params.tx_gain)

                try
                    % Update progress
                    progress_it = (y - 1)*length(real_params.chip_no)*length(real_params.tx_gain) + ...
                        (z - 1)*length(real_params.tx_gain) + x;
                
                    % Reset time variable so Bob knows when Alice finishes transmitting
                    % Pause above give's Bob enough time see that Alice as finished
                    current_time_tx_alice = 0;
                    save('alice_time.mat','current_time_tx_alice')
                    
                    % Wait for Bob to stop listening to noise
                    % Noise function shouldn't throw an error
                    disp('Pause to wait for Bob to finish listening to noise.')
                    current_time_noise_bob = 0;
                    for count = 1 : real_params.count_max
                       if current_time_noise_bob ~= 0
                         disp(['Noise reception for Transmission ', num2str(y), '.', num2str(z), '.', num2str(x), ' has ended'])
                         break; % terminate execution of this smaller for-loop
                       end
                       pause(0.25) % Can load pretty often thanks to try/catch
                       try
                           load('bob_time.mat','current_time_noise_bob')
                       catch
                           continue
                       end
                    end
                    
                    % Start timer for one transmission
                    t_start_tx(y,x,z) = tic;
                    
                    % Display message in command window
                    disp(['Transmission ', num2str(y), '.', num2str(z), '.', num2str(x), ' has started'])
                    % Transmit function
                    [txSig(y,x,z), current_time_tx_alice] = tx_alice(x, z, frame_size(1,z), real_params);
                
                    % Non-zero current time
                    if current_time_tx_alice ~= 0
                        % Display message in command window
                        disp(['Transmission ', num2str(y), '.', num2str(z), '.', num2str(x), ' has ended'])
                        % Save current time so Bob knows when Alice is finished
                        save('alice_time.mat','current_time_tx_alice')
                    end   
                    
                    % Record time taken for one transmission
                    t_end_tx(y,x,z) = toc(t_start_tx(y,x,z));
                
                    % Pause so Bob has time to realise Alice's transmission has ended
                    % About to wait for noise reception anyway
                    pause(1)

                catch ME % We have encountered an error above

                    % Let's Bob know that Alice has finished (non-zero time)
                    current_time_tx_alice = 1;
                    save('alice_time.mat','current_time_tx_alice')

                    % Register the fact we encountered an error
                    alice_err_log(progress_it).err_msg = ME.message;
                    alice_err_log(progress_it).err_func1 = ME.stack(2).name;
                    alice_err_log(progress_it).err_func2 = ME.stack(1).name;
                    % Simplest way to access precise iteration
                    alice_err_log(progress_it).err_y = y;
                    alice_err_log(progress_it).err_x = x;
                    alice_err_log(progress_it).err_z = z;

                    % Wait for Bob to finish listening
                    disp('Error encountered.')
                    disp('Pause to wait for Bob to finish listening')
                    current_time_rx_bob = 0;
                    for count = 1 : real_params.count_max
                       if current_time_rx_bob ~= 0
                         disp(['Reception of Transmission ', num2str(y), '.', num2str(z), '.', num2str(x),' has ended'])
                         break; % terminate execution of this smaller for-loop
                       end
                       pause(0.25) % Can load pretty often thanks to try/catch
                       try
                           load('bob_time.mat','current_time_rx_bob')
                       catch
                           continue
                       end
                    end

                    % Pause so Bob has time to realise Alice has failed
                    % About to wait for noise reception anyway
                    pause(1)

                end
        
            end
    
        end
    
    end

    %% ERROR HANDLING

    % Save Alice's error results so Bob is aware
    save('alice_errors.mat', 'alice_err_log')

    %% SAVE

    % Save results to be safe, use -v7.3 MAT-file version to store large data
    save('alice_results_.mat', ...
        'txSig','t_end_tx', ...
        'alice_err_log', ...
        'real_params',...
        '-v7.3')
    disp('**WARNING**: DO NOT RUN PROGRAM AGAIN WITHOUT CHANGING RESULTS FILE NAME!')

else
    % End program if less than two radios are connected
    disp('Less than two radios are connected')
    return
end