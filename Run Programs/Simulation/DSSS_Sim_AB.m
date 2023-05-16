% Program to simulate transmission between Alice and Bob
% With Willie listening

clear

%% INITIALISE
% Transmitter parameter structure
sim_params = sim_init;
% Preallocate output vectors
% Using NaN vectors helps if an experiment is cut short
% Signal vectors
signal_template = cell(sim_params.data_reps, length(sim_params.SNR_dB), length(sim_params.chip_no));
frame_size = NaN(1, length(sim_params.chip_no));
% Signal powers
results_template = NaN(sim_params.data_reps, length(sim_params.SNR_dB), length(sim_params.chip_no));
noise_pwr = results_template;
txSig_pwr = results_template;
% SNR & BER Results
SNR_dB = repmat(sim_params.SNR_dB, sim_params.data_reps, 1, length(sim_params.chip_no));
BER_theo = results_template;
rxSig_BER = repmat({NaN(3,1)},sim_params.data_reps, length(sim_params.SNR_dB), length(sim_params.chip_no));
% Willie's hypothesis tests
LLR_noise = results_template;
LLR_sig = results_template;

% Start timer for whole program
t_start_prog = tic;

% Vectors and cell arrays that can only be defined with (varying) frame size
for z = 1 : length(sim_params.chip_no)

    % Spreading will dictate frame size (and thus frame time)
    % Can be an even positive integer from 2 to 16,777,216
    % Using values less than 3660 can yield poor performance
    frame_size(1,z) = (sim_params.header_len + sim_params.payload_len)*sim_params.chip_no(z)/sim_params.bps;

    % Preallocate 3D cell array of signal vectors
    % frame_size changes with each "page"
    signal_template(:,:,z) = repmat({complex(zeros(frame_size(1,z) * sim_params.interpol, 1), 0)}, ...
        sim_params.data_reps, length(sim_params.SNR_dB));
    noise_sig = signal_template;
    rx_sig = signal_template;

end

% Waitbar to view progress of program
w_bar = waitbar(0, 'First transmission');
program_it = length(sim_params.SNR_dB) * sim_params.data_reps * length(sim_params.chip_no);
wbar_text = ['Transmission No %d.%d.%d out of %d.%d.%d completed\n' ...
            'Progress: %d %%'];

%% ITERATE
% Iterate over different chip and receiver gain values, repeatedly

for y = 1 : sim_params.data_reps

    for z = 1 : length(sim_params.chip_no)

        for x = 1 : length(sim_params.SNR_dB)
        
            [rx_sig(y,x,z), noise_sig(y,x,z), ...
                SNR_dB(y,x,z), txSig_pwr(y,x,z), noise_pwr(y,x,z), ...
                BER_theo(y,x,z), rxSig_BER{y,x,z}] = tx_rx_sim(x, z, sim_params);
    
            % Willie calculates LLR value for each signal
            LLR_noise(y,x,z) = LLR_calc(noise_sig{y,x,z}, noise_pwr(y,x,z), txSig_pwr(y,x,z), sim_params);
            LLR_sig(y,x,z) = LLR_calc(rx_sig{y,x,z}, noise_pwr(y,x,z), txSig_pwr(y,x,z), sim_params);
    
            % Update waitbar
            progress_it = (y - 1)*length(sim_params.chip_no)*length(sim_params.SNR_dB) + (z - 1)*length(sim_params.SNR_dB) + x;
            progress_perc = floor(progress_it/program_it*100);
            w_bar = waitbar(progress_it/program_it, w_bar, ...
                sprintf(wbar_text,y,z,x, ...
                sim_params.data_reps, ...
                length(sim_params.chip_no), ...
                length(sim_params.SNR_dB), ...
                progress_perc));
            pause(0.1);
    
        end

    end

end

% Record time taken for whole program
t_end_prog = toc(t_start_prog);
avrg_it_time = t_end_prog/program_it;

% Save results to be safe, use -v7.3 MAT-file version to store large data
save('sim_results_.mat', ...
    'txSig_pwr','SNR_dB', ...
    'rx_sig','rxSig_BER','BER_theo', ...
    'noise_sig','noise_pwr', ...
    'LLR_noise','LLR_sig',...
    'avrg_it_time', ...
    'frame_size', ...
    'sim_params',...
    '-v7.3')
disp('**WARNING**: DO NOT RUN PROGRAM AGAIN WITHOUT CHANGING RESULTS FILE NAME!')
