% Program to receive Alice's transmission

% Make sure you've opened MATLAB with administrator privileges

clear

% Store and check information regarding attached radios
connected_radios_rx = findPlutoRadio;
if length(connected_radios_rx) == 3
    %% INITIALISE
    % Transmitter parameter structure
    sim_params = sim_init;
    % Preallocate output vectors
    % Using NaN vectors helps if an experiment is cut short
    % Signal vectors
    signal_template = cell(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no));
    frame_size = NaN(1, length(sim_params.chip_no));
    % Signal powers
    results_template = NaN(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no));
    noise_pwr = results_template;
    txSig_pwr = results_template;
    % SNR & BER Results
    SNR_dB = results_template;
    SNR_dB_theo = results_template;
    BER_theo = results_template;
    rxSig_BER = repmat({NaN(3,1)},sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no));
    % Noise characteristics
    %mu_sig = repmat({NaN(1,2)},sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no));
    %H0_lil = results_template;
    noise_fig = results_template;
    noise_flr = results_template;
    % Frequency correction factor
    freq_corr = NaN(sim_params.data_reps, length(sim_params.chip_no));
    
    % Preallocate stopwatch vectors
    time_template = zeros(sim_params.data_reps, length(sim_params.tx_gain), length(sim_params.chip_no), 'uint64');
    t_start_noise = time_template;
    t_end_noise = time_template;
    t_start_rx = time_template;
    t_end_rx = time_template;
    
    %% SYNCHRONISE
    
    % Reset stopwatch variable
    % Helps if previous run was cut short by error
    current_time_noise_bob = 0;
    save('noise_time_bob.mat','current_time_noise_bob')
    
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
    
    % Waitbar to view progress of program
    w_bar = waitbar(0, 'First transmission');
    program_it = length(sim_params.tx_gain) * sim_params.data_reps * length(sim_params.chip_no);
    wbar_text = ['Transmission No %d.%d.%d out of %d.%d.%d completed\n' ...
                'Progress: %d %%'];
    
    %% ITERATE
    % Iterate over different chip and receiver gain values, repeatedly
    
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
            noise_sig = signal_template;
            rx_sig = signal_template;
    
            % Determine frequency correction factor for every new chip value
            freq_corr(y,z) = freq_correct(frame_size(1,z), sim_params);
    
            for x = 1 : length(sim_params.tx_gain)
        
                % Reset for syncing up transmission/reception inside each loop
                % Need to reset current_time_noise_bob variable to zero for Alice
                current_time_noise_bob = 0;
                save('noise_time_bob.mat','current_time_noise_bob')
            
                %% NOISE RECEPTION
            
                % Start timer for one noise loop
                t_start_noise(y,x,z) = tic;
            
                disp(['Obtaining noise sample for Transmission ', num2str(y), '.', num2str(z), '.', num2str(x)])
                [noise_sig(y,x,z), noise_pwr(y,x,z), SNR_dB_theo(y,x,z), noise_fig(y,x,z), noise_flr(y,x,z), current_time_noise_bob] = ...
                    rx_noise(x, frame_size(1,z), freq_corr(y,z), sim_params.Bob_Noise_ID, sim_params);
        
                if current_time_noise_bob ~= 0
                    disp(['Noise reception for Transmission ', num2str(y), '.', num2str(z), '.', num2str(x),' has ended'])
                    % Save current time so Alice knows when Bob is finished
                    save('noise_time_bob.mat','current_time_noise_bob')
                end 
            
                % Record time taken for one loop
                t_end_noise(y,x,z) = toc(t_start_noise(y,x,z));
            
                % Pause so Alice has time to realise Bob has finished receiving noise
                pause(0.75) %1.5
            
                %% SIGNAL RECEPTION
                
                disp('Pause to ensure Alice is already transmitting.')
                pause(0.75) %1.5
            
                % Start timer for one reception loop
                t_start_rx(y,x,z) = tic;
            
                disp(['Listening for Transmission ', num2str(y), '.', num2str(z), '.', num2str(x)])
                [SNR_dB(y,x,z), BER_theo(y,x,z), rxSig_BER{y,x,z}, rx_sig(y,x,z), txSig_pwr(y,x,z), current_time_rx_bob] = ...
                rx_bob(z, frame_size(1,z),  noise_pwr(y,x,z), freq_corr(y,z), sim_params);
        
                if current_time_rx_bob ~= 0
                    disp(['Reception of Transmission ', num2str(y), '.', num2str(z), '.', num2str(x),' has ended'])
                end   
            
                % Record time taken for one reception loop
                t_end_rx(y,x,z) = toc(t_start_rx(y,x,z));
        
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
        
                % Update waitbar
                progress_it = (y - 1)*length(sim_params.tx_gain) + (z - 1)*length(sim_params.chip_no) + x;
                progress_perc = floor(progress_it/program_it*100);
                w_bar = waitbar(progress_it/program_it, w_bar, ...
                    sprintf(wbar_text,y,z,x, ...
                    sim_params.data_reps, ...
                    length(sim_params.chip_no), ...
                    length(sim_params.tx_gain), ...
                    progress_perc));
                pause(0.1);
        
            end
    
        end
    
    end
    
    % Save results to be safe, use -v7.3 MAT-file version to store large data
    save('bob_results_.mat', ...
        'txSig_pwr','SNR_dB','SNR_dB_theo', ...
        'rx_sig','rxSig_BER','BER_theo', ...
        'noise_sig','noise_pwr','noise_fig','noise_flr', ...
        't_end_rx','t_end_noise', ...
        'freq_corr', ...
        'frame_size', ...
        'sim_params',...
        '-v7.3')
    disp('**WARNING**: DO NOT RUN PROGRAM AGAIN WIHTOUT CHANGING RESULTS FILE NAME!')

else
    % End program if less than two radios are connected
    disp('Less than three radios are connected')
    return
end