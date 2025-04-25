% MIMO System with Enhanced DOA Estimation, Improved Doppler Compensation and Artificial Noise
% This simulation evaluates BER performance under various conditions with security features
close all;
clc;

%% Simulation Parameters
SNR_dB = 0:5:30;         % SNR range in dB
SNR_linear = 10.^(SNR_dB/10); % SNR in linear scale
nBits = 5000;            % Number of bits per SNR point
nBlocks = 1000;          % Number of channel blocks
modulation_order = 4;    % QPSK
doppler_fd_levels = [0, 0.001, 0.003, 0.01];  % Normalized Doppler frequencies
Nt = 2;                  % Number of transmit antennas
Nr = 2;                  % Number of receive antennas
Ns = 1;                  % Number of data streams
AN_power_ratio = 0.5;    % Power ratio allocated to AN

%% Channel and Signal Parameters
carrier_freq = 2.4e9;    % 2.4 GHz
sample_time = 1e-6;      % 1 μs
doppler_phase_impact = 2*pi;  % Base impact of Doppler on phase

% Channel Coding Parameters
trellis = poly2trellis(9, [561 753]);  % Rate 1/2 convolutional code
code_rate = 1/2;                       % Code rate
tracebackLength = 60;                  % Traceback length
termination_length = 8;                % Termination zeros

% Pilot structure parameters
base_pilot_spacing = 20;   % Base pilot spacing (for fd = 0)
min_pilot_spacing = 5;     % Minimum pilot spacing
pilot_values = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2); % Normalized QPSK constellation

% ULA Parameters for DOA Estimation
antenna_spacing = 0.5;     % Spacing in wavelengths
n_angles = 181;            % Number of angles to scan
angle_range = linspace(0, 180, n_angles); % Angles in degrees
angle_range_rad = angle_range * pi/180;

% Eavesdropper parameters
eavesdropper_angles = [30, 150]; % Degrees

%% Preallocate arrays for results
BER_no_comp = zeros(length(doppler_fd_levels), length(SNR_dB));
BER_with_comp = zeros(length(doppler_fd_levels), length(SNR_dB));
BER_with_kalman = zeros(length(doppler_fd_levels), length(SNR_dB));
BER_eve = zeros(length(doppler_fd_levels), length(SNR_dB));

actual_BER_no_comp = zeros(length(doppler_fd_levels), length(SNR_dB));
actual_BER_with_comp = zeros(length(doppler_fd_levels), length(SNR_dB));
actual_BER_with_kalman = zeros(length(doppler_fd_levels), length(SNR_dB));
actual_BER_eve = zeros(length(doppler_fd_levels), length(SNR_dB));

doppler_est_values = zeros(length(doppler_fd_levels), 1);
doa_estimates = zeros(length(doppler_fd_levels), length(SNR_dB), nBlocks, 2);
doa_errors = zeros(length(doppler_fd_levels), length(SNR_dB));

phase_mse_basic = zeros(length(doppler_fd_levels), length(SNR_dB));
phase_mse_kalman = zeros(length(doppler_fd_levels), length(SNR_dB));

true_angles = [30, 120]; % Two signal sources

fprintf('Testing MIMO System with Artificial Noise, Enhanced Doppler Compensation, and DOA Estimation\n');
fprintf('-----------------------------------------------------------------------------------\n');

%% Main simulation loop
% NEW: Ensure doppler_fd_levels is defined
if ~exist('doppler_fd_levels', 'var')
    doppler_fd_levels = [0, 0.001, 0.003, 0.01];
    warning('doppler_fd_levels was undefined. Reinitialized to default values.');
end

for fd_idx = 1:length(doppler_fd_levels)
    doppler_fd = doppler_fd_levels(fd_idx);
    
    fprintf('Processing Doppler frequency fd = %g\n', doppler_fd);
    
    doppler_est_values(fd_idx) = 0;
    
    % Exponential decay for pilot spacing
    if doppler_fd == 0
        pilot_spacing = base_pilot_spacing;
    else
        pilot_spacing = round(base_pilot_spacing - (base_pilot_spacing - min_pilot_spacing) * (1 - exp(-doppler_fd/0.005)));
    end
    
    fprintf('Using pilot spacing: %d symbols\n', pilot_spacing);

    for snr_idx = 1:length(SNR_dB)
        current_snr_db = SNR_dB(snr_idx);
        current_snr = SNR_linear(snr_idx);
        noisePower = 1/current_snr;
        
        fprintf('  SNR = %d dB: ', current_snr_db);
        
        if current_snr_db > 10
            effective_pilot_spacing = min(pilot_spacing * 2, 30);
        else
            effective_pilot_spacing = pilot_spacing;
        end
        
        errors_no_comp = 0;
        errors_with_comp = 0;
        errors_with_kalman = 0;
        errors_eve = 0;
        total_bits_no_comp = 0;
        total_bits_with_comp = 0;
        total_bits_with_kalman = 0;
        total_bits_eve = 0;
        
        phase_errors_basic = 0;
        phase_errors_kalman = 0;
        phase_samples = 0;
        
        for blk = 1:nBlocks
            txDataBits = randi([0 1], floor((nBits * code_rate) - termination_length), 1);
            txBits = [txDataBits; zeros(termination_length, 1)];
            codedBits = convenc(txBits, trellis);
            intrlvrIndices = randperm(length(codedBits));
            intlvdBits = codedBits(intrlvrIndices);
            txSymbols = qammod(intlvdBits, modulation_order, 'InputType', 'bit', 'UnitAveragePower', true);
            nSymbolsPerStream = floor(length(txSymbols)/Ns);
            txSymbolsReshaped = reshape(txSymbols(1:nSymbolsPerStream*Ns), nSymbolsPerStream, Ns);
            
            nDataSymbols = nSymbolsPerStream;
            % NEW: Ensure at least 2 pilots
            nPilots = max(2, floor(nDataSymbols / effective_pilot_spacing));
            nTotalSymbols = nDataSymbols + nPilots;
            txSymbolsWithPilots = zeros(nTotalSymbols, Ns);
            pilot_positions = [];
            pilot_values_used = [];
            data_idx = 1;
            pilot_idx = 1;
            
            for i = 1:nTotalSymbols
                if mod(i-1, effective_pilot_spacing+1) == 0 && pilot_idx <= nPilots
                    current_pilot = pilot_values(mod(pilot_idx-1, length(pilot_values))+1);
                    txSymbolsWithPilots(i, :) = current_pilot;
                    pilot_positions = [pilot_positions, i];
                    pilot_values_used = [pilot_values_used, current_pilot];
                    pilot_idx = pilot_idx + 1;
                else
                    if data_idx <= nDataSymbols
                        txSymbolsWithPilots(i, :) = txSymbolsReshaped(data_idx, :);
                        data_idx = data_idx + 1;
                    end
                end
            end

            signal_angles_rad = true_angles * pi/180;
            steering_vectors = zeros(Nr, length(signal_angles_rad));
            for i = 1:length(signal_angles_rad)
                for n = 0:Nr-1
                    steering_vectors(n+1, i) = exp(-1j * 2 * pi * n * antenna_spacing * cos(signal_angles_rad(i)));
                end
            end
            H_base_coeffs = (randn(length(signal_angles_rad), Nt) + 1j*randn(length(signal_angles_rad), Nt))/sqrt(2);
            H_base = steering_vectors * H_base_coeffs;
            [U_base, S_base, V_base] = svd(H_base);
            W = V_base(:, 1:Ns);
            G = U_base(:, 1:Ns)';
            
            num_eve = length(eavesdropper_angles);
            He = zeros(num_eve, Nt);
            for e = 1:num_eve
                angle_rad = eavesdropper_angles(e) * pi/180;
                for tx = 1:Nt
                    He(e, tx) = exp(-1j * 2 * pi * antenna_spacing * (tx-1) * sin(angle_rad));
                end
            end
            He = He ./ sqrt(sum(abs(He).^2, 2)) .* (randn(num_eve, Nt) + 1i*randn(num_eve, Nt))/sqrt(2);
            
            [U, S, V] = svd(H_base);
            signal_subspace = V(:, 1:Ns);
            AN_subspace = V(:, Ns+1:end);
            if isempty(AN_subspace)
                AN = (randn(Nt, size(txSymbolsWithPilots, 1)) + 1i*randn(Nt, size(txSymbolsWithPilots, 1)))/sqrt(2);
                proj = (signal_subspace * signal_subspace') * AN;
                AN = AN - proj;
                AN = AN * sqrt(AN_power_ratio / mean(sum(abs(AN).^2, 1)));
            else
                AN_coeffs = (randn(size(AN_subspace, 2), size(txSymbolsWithPilots, 1)) + ...
                          1i*randn(size(AN_subspace, 2), size(txSymbolsWithPilots, 1)))/sqrt(2);
                AN = AN_subspace * AN_coeffs;
                AN = AN * sqrt(AN_power_ratio / mean(sum(abs(AN).^2, 1)));
            end
            
            tx_signal_precoded = W * txSymbolsWithPilots.';
            signal_power = 1 - AN_power_ratio;
            tx_signal_precoded = sqrt(signal_power) * tx_signal_precoded;
            tx_signal_with_AN = tx_signal_precoded + AN;
            
            rx_signal_no_comp = zeros(Nr, size(tx_signal_with_AN, 2));
            rx_signal_eve = zeros(num_eve, size(tx_signal_with_AN, 2));
            
            for t = 1:size(tx_signal_with_AN, 2)
                if doppler_fd > 0
                    doppler_phase = 2*pi*doppler_fd*t;
                    H_t = H_base * exp(1i * doppler_phase);
                else
                    H_t = H_base;
                end
                rx_t = H_t * tx_signal_with_AN(:, t);
                rx_signal_no_comp(:, t) = rx_t + sqrt(noisePower/2) * (randn(Nr, 1) + 1i*randn(Nr, 1));
                for e = 1:num_eve
                    rx_eve = He(e,:) * tx_signal_with_AN(:, t);
                    rx_signal_eve(e, t) = rx_eve + sqrt(noisePower/2) * (randn + 1i*randn);
                end
            end
            rx_processed_raw = G * rx_signal_no_comp;

            %% Enhanced Direction of Arrival (DOA) Estimation
            if mod(blk, 10) == 0 || blk == 1
                current_snr_db = SNR_dB(snr_idx);
                [est_angles, P_MUSIC] = improved_doa_estimation(rx_signal_no_comp, antenna_spacing, n_angles, angle_range_rad, Ns);
                doa_estimates(fd_idx, snr_idx, blk, :) = est_angles;
                
                if est_angles(1) > 0 && est_angles(2) > 0
                    error1a = abs(est_angles(1) - true_angles(1)) + abs(est_angles(2) - true_angles(2));
                    error1b = abs(est_angles(1) - true_angles(2)) + abs(est_angles(2) - true_angles(1));
                    if error1a < error1b
                        cur_error = error1a / 2;
                    else
                        cur_error = error1b / 2;
                    end
                elseif est_angles(1) > 0
                    cur_error = min(abs(est_angles(1) - true_angles));
                else
                    cur_error = 90;
                end
                if blk == 1
                    doa_errors(fd_idx, snr_idx) = cur_error;
                else
                    doa_errors(fd_idx, snr_idx) = (doa_errors(fd_idx, snr_idx) * (blk/10-1) + cur_error) / (blk/10);
                end
                if blk == 1 && snr_idx == length(SNR_dB)
                    fprintf('\nDOA Estimation Results (fd = %.4f, SNR = %d dB):\n', doppler_fd, current_snr_db);
                    fprintf('  True angles: [%s]\n', num2str(true_angles));
                    fprintf('  Estimated angles: [%s]\n', num2str(est_angles));
                    if est_angles(1) > 0 && est_angles(2) > 0
                        if error1a < error1b
                            fprintf('  Mean angular error: %.2f degrees\n\n', error1a/2);
                        else
                            fprintf('  Mean angular error: %.2f degrees\n\n', error1b/2);
                        end
                    elseif est_angles(1) > 0
                        fprintf('  Angular error: %.2f degrees (only one angle detected)\n\n', min(abs(est_angles(1) - true_angles)));
                    else
                        fprintf('  No angles detected\n\n');
                    end
                end
            end
            
            %% Doppler Shift Estimation (Legitimate Receiver)
            estimated_doppler_fd = 0;
            if ~isempty(pilot_positions) && length(pilot_positions) >= 2
                pilot_phases = zeros(length(pilot_positions), Ns);
                for i = 1:length(pilot_positions)
                    pos = pilot_positions(i);
                    actual_pilot = pilot_values_used(i);
                    for n = 1:Ns
                        pilot_phases(i, n) = angle(rx_processed_raw(n, pos) / actual_pilot);
                    end
                end
                for n = 1:Ns
                    pilot_phases(:, n) = unwrap(pilot_phases(:, n));
                end
                phase_diffs = zeros(length(pilot_positions)-1, Ns);
                time_diffs = zeros(length(pilot_positions)-1, 1);
                for i = 1:length(pilot_positions)-1
                    time_diffs(i) = pilot_positions(i+1) - pilot_positions(i);
                    for n = 1:Ns
                        phase_diffs(i, n) = pilot_phases(i+1, n) - pilot_phases(i, n);
                    end
                end
                fd_estimates = phase_diffs ./ repmat(time_diffs, 1, Ns) / (2*pi);
                estimated_doppler_fd = median(fd_estimates(:));
                if abs(estimated_doppler_fd) > 1.5 * max(doppler_fd_levels)
                    estimated_doppler_fd = sign(estimated_doppler_fd) * max(doppler_fd_levels);
                end
            end
            if snr_idx == 1 && blk == 1
                fprintf('True Doppler fd = %.6f, Estimated fd = %.6f\n', doppler_fd, estimated_doppler_fd);
            end
            if SNR_dB(snr_idx) >= 15 && SNR_dB(snr_idx) <= 20
                doppler_est_values(fd_idx) = doppler_est_values(fd_idx) + estimated_doppler_fd / nBlocks / length(find(SNR_dB >= 15 & SNR_dB <= 20));
            end

            %% 1. Processing WITHOUT Doppler compensation
            rx_processed_no_comp = G * rx_signal_no_comp;
            rx_data_no_comp = zeros(nDataSymbols, Ns);
            data_idx = 1;
            for i = 1:size(rx_processed_no_comp, 2)
                if ~ismember(i, pilot_positions) && data_idx <= nDataSymbols
                    rx_data_no_comp(data_idx, :) = rx_processed_no_comp(:, i).';
                    data_idx = data_idx + 1;
                end
            end
            rxCodedBits_no_comp = qamdemod(rx_data_no_comp(:), modulation_order, 'OutputType', 'bit', 'UnitAveragePower', true);
            deintlvdBits_no_comp = zeros(size(rxCodedBits_no_comp));
            deintlvdBits_no_comp(intrlvrIndices) = rxCodedBits_no_comp;
            rxBits_no_comp = vitdec(deintlvdBits_no_comp, trellis, tracebackLength, 'term', 'hard');
            rxDataBits_no_comp = rxBits_no_comp(1:end-termination_length);
            num_errors_no_comp = sum(rxDataBits_no_comp ~= txDataBits(1:min(length(rxDataBits_no_comp), length(txDataBits))));
            errors_no_comp = errors_no_comp + num_errors_no_comp;
            total_bits_no_comp = total_bits_no_comp + min(length(rxDataBits_no_comp), length(txDataBits));
            
            %% 2. Processing WITH Basic Doppler compensation
            if doppler_fd < 0.0005
                rx_processed_comp = rx_processed_raw;
            else
                rx_processed_comp = zeros(size(rx_processed_raw));
                if ~isempty(pilot_positions) && length(pilot_positions) >= 2
                    rx_pilots = zeros(length(pilot_positions), Ns);
                    rx_pilots_actual = zeros(length(pilot_positions), Ns);
                    for i = 1:length(pilot_positions)
                        pos = pilot_positions(i);
                        rx_pilots(i, :) = rx_processed_raw(:, pos).';
                        rx_pilots_actual(i, :) = pilot_values_used(i);
                    end
                    pilot_phases = zeros(size(rx_pilots));
                    for i = 1:size(rx_pilots, 1)
                        for j = 1:Ns
                            pilot_phases(i, j) = angle(rx_pilots(i, j) / rx_pilots_actual(i, j));
                        end
                    end
                    all_phases = zeros(size(rx_processed_raw, 2), Ns);
                    for n = 1:Ns
                        unwrapped_phases = unwrap(pilot_phases(:, n));
                        if length(unwrapped_phases) >= 3
                            unwrapped_phases = movmean(unwrapped_phases, 3);
                        end
                        all_phases(:, n) = interp1(pilot_positions, unwrapped_phases, 1:size(rx_processed_raw, 2), 'pchip', 'extrap');
                    end
                    for t = 1:size(rx_processed_raw, 2)
                        for n = 1:Ns
                            rx_processed_comp(n, t) = rx_processed_raw(n, t) * exp(-1i * all_phases(t, n));
                        end
                    end
                    % Compute phase MSE for basic compensation
                    if doppler_fd > 0
                        true_phases = 2*pi*doppler_fd*(1:size(rx_processed_raw, 2));
                        for n = 1:Ns
                            phase_errors_basic = phase_errors_basic + sum((all_phases(:, n)' - true_phases).^2);
                        end
                        phase_samples = phase_samples + size(rx_processed_raw, 2) * Ns;
                    end
                else
                    rx_processed_comp = rx_processed_raw;
                end
            end
            
            rx_data_comp = zeros(nDataSymbols, Ns);
            data_idx = 1;
            for i = 1:size(rx_processed_comp, 2)
                if ~ismember(i, pilot_positions) && data_idx <= nDataSymbols
                    rx_data_comp(data_idx, :) = rx_processed_comp(:, i).';
                    data_idx = data_idx + 1;
                end
            end
            rxCodedBits_with_comp = qamdemod(rx_data_comp(:), modulation_order, 'OutputType', 'bit', 'UnitAveragePower', true);
            deintlvdBits_with_comp = zeros(size(rxCodedBits_with_comp));
            deintlvdBits_with_comp(intrlvrIndices) = rxCodedBits_with_comp;
            rxBits_with_comp = vitdec(deintlvdBits_with_comp, trellis, tracebackLength, 'term', 'hard');
            rxDataBits_with_comp = rxBits_with_comp(1:end-termination_length);
            num_errors_with_comp = sum(rxDataBits_with_comp ~= txDataBits(1:min(length(rxDataBits_with_comp), length(txDataBits))));
            errors_with_comp = errors_with_comp + num_errors_with_comp;
            total_bits_with_comp = total_bits_with_comp + min(length(rxDataBits_with_comp), length(txDataBits));

            %% 3. Processing WITH Enhanced Kalman filtering
            if doppler_fd < 0.0005
                rx_processed_kalman = rx_processed_raw;
            else
                rx_processed_kalman = zeros(size(rx_processed_raw));
                if ~isempty(pilot_positions) && length(pilot_positions) >= 2
                    % Enhanced Kalman parameters
                    x_est = [0; 2*pi*estimated_doppler_fd; 0];
                    P_est = [0.01, 0, 0; 0, 0.001, 0; 0, 0, 0.0001];
                    Q = [5e-4, 0, 0; 0, 5e-5, 0; 0, 0, 5e-6] * sqrt(max(doppler_fd, 0.0001) / 0.01);
                    R = min(0.02 * (1 + 1/max(current_snr, 0.1)), 0.1);
                    dt = 1;
                    F = [1, dt, 0.5*dt^2; 0, 1, dt; 0, 0, 1];
                    H = [1, 0, 0];
                    
                    % Least-squares initialization
                    if length(pilot_positions) >= 3
                        pilot_phases = zeros(length(pilot_positions), 1);
                        for i = 1:length(pilot_positions)
                            pos = pilot_positions(i);
                            pilot_phases(i) = angle(rx_processed_raw(1, pos) / pilot_values_used(i));
                        end
                        pilot_phases = unwrap(pilot_phases);
                        t_pilot = pilot_positions';
                        t_mean = mean(t_pilot);
                        phi_mean = mean(pilot_phases);
                        v_est = sum((t_pilot - t_mean) .* (pilot_phases - phi_mean)) / sum((t_pilot - t_mean).^2);
                        phi_est = pilot_phases(1);
                        x_est = [phi_est; v_est; 0];
                    end
                    
                    kalman_phases = zeros(size(rx_processed_raw, 2), Ns);
                    for n = 1:Ns
                        x_curr = x_est;
                        P_curr = P_est;
                        all_measurements = zeros(length(pilot_positions), 1);
                        for i = 1:length(pilot_positions)
                            pos = pilot_positions(i);
                            actual_pilot = pilot_values_used(i);
                            all_measurements(i) = angle(rx_processed_raw(n, pos) / actual_pilot);
                        end
                        all_measurements = unwrap(all_measurements);
                        
                        for t = 1:size(rx_processed_raw, 2)
                            x_pred = F * x_curr;
                            P_pred = F * P_curr * F' + Q;
                            if ismember(t, pilot_positions)
                                pilot_idx = find(pilot_positions == t);
                                measured_phase = all_measurements(pilot_idx);
                                S = H * P_pred * H' + R;
                                K = P_pred * H' / S;
                                innovation = measured_phase - H * x_pred;
                                x_curr = x_pred + K * innovation;
                                I_KH = eye(size(P_curr)) - K * H;
                                P_curr = I_KH * P_pred * I_KH' + K * R * K';
                                % NEW: Enforce positive definiteness
                                [V, D] = eig(P_curr);
                                D = max(D, 1e-6);
                                P_curr = V * D * V';
                            else
                                x_curr = x_pred;
                                P_curr = P_pred;
                            end
                            kalman_phases(t, n) = x_curr(1);
                        end
                    end
                    for t = 1:size(rx_processed_raw, 2)
                        for n = 1:Ns
                            rx_processed_kalman(n, t) = rx_processed_raw(n, t) * exp(-1i * kalman_phases(t, n));
                        end
                    end
                    % Compute phase MSE for Kalman
                    if doppler_fd > 0
                        true_phases = 2*pi*doppler_fd*(1:size(rx_processed_raw, 2));
                        for n = 1:Ns
                            phase_errors_kalman = phase_errors_kalman + sum((kalman_phases(:, n)' - true_phases).^2);
                        end
                        phase_samples = phase_samples + size(rx_processed_raw, 2) * Ns;
                    end
                else
                    rx_processed_kalman = rx_processed_raw;
                end
            end
            
            rx_data_kalman = zeros(nDataSymbols, Ns);
            data_idx = 1;
            for i = 1:size(rx_processed_kalman, 2)
                if ~ismember(i, pilot_positions) && data_idx <= nDataSymbols
                    rx_data_kalman(data_idx, :) = rx_processed_kalman(:, i).';
                    data_idx = data_idx + 1;
                end
            end
            rxCodedBits_with_kalman = qamdemod(rx_data_kalman(:), modulation_order, 'OutputType', 'bit', 'UnitAveragePower', true);
            deintlvdBits_with_kalman = zeros(size(rxCodedBits_with_kalman));
            deintlvdBits_with_kalman(intrlvrIndices) = rxCodedBits_with_kalman;
            rxBits_with_kalman = vitdec(deintlvdBits_with_kalman, trellis, tracebackLength, 'term', 'hard');
            rxDataBits_with_kalman = rxBits_with_kalman(1:end-termination_length);
            num_errors_kalman = sum(rxDataBits_with_kalman ~= txDataBits(1:min(length(rxDataBits_with_kalman), length(txDataBits))));
            errors_with_kalman = errors_with_kalman + num_errors_kalman;
            total_bits_with_kalman = total_bits_with_kalman + min(length(rxDataBits_with_kalman), length(txDataBits));

            %% Process Eavesdropper Signal
            % NEW: Doppler estimation for eavesdropper
            estimated_doppler_fd_eve = 0;
            if doppler_fd > 0 && ~isempty(pilot_positions) && length(pilot_positions) >= 2
                pilot_phases_eve = zeros(length(pilot_positions), 1);
                for i = 1:length(pilot_positions)
                    pos = pilot_positions(i);
                    actual_pilot = pilot_values_used(i);
                    pilot_phases_eve(i) = angle(rx_signal_eve(1, pos) / actual_pilot);
                end
                pilot_phases_eve = unwrap(pilot_phases_eve);
                phase_diffs_eve = diff(pilot_phases_eve);
                time_diffs_eve = diff(pilot_positions');
                fd_estimates_eve = phase_diffs_eve ./ time_diffs_eve / (2*pi);
                estimated_doppler_fd_eve = median(fd_estimates_eve);
                if abs(estimated_doppler_fd_eve) > 1.5 * max(doppler_fd_levels)
                    estimated_doppler_fd_eve = sign(estimated_doppler_fd_eve) * max(doppler_fd_levels);
                end
            end
            
            rx_eve = rx_signal_eve(1, :);
            rx_processed_eve = zeros(1, length(rx_eve));
            for t = 1:length(rx_eve)
                doppler_phase = 2*pi*estimated_doppler_fd_eve*t;
                rx_processed_eve(t) = rx_eve(t) * exp(-1i * doppler_phase);
            end
            rx_data_eve = zeros(nDataSymbols, 1);
            data_idx = 1;
            for i = 1:length(rx_processed_eve)
                if ~ismember(i, pilot_positions) && data_idx <= nDataSymbols
                    rx_data_eve(data_idx) = rx_processed_eve(i);
                    data_idx = data_idx + 1;
                end
            end
            rxCodedBits_eve = qamdemod(rx_data_eve, modulation_order, 'OutputType', 'bit', 'UnitAveragePower', true);
            deintlvdBits_eve = zeros(size(rxCodedBits_eve));
            deintlvdBits_eve(intrlvrIndices(1:length(rxCodedBits_eve))) = rxCodedBits_eve;
            rxBits_eve = vitdec(deintlvdBits_eve, trellis, tracebackLength, 'term', 'hard');
            rxDataBits_eve = rxBits_eve(1:end-termination_length);
            num_errors_eve = sum(rxDataBits_eve ~= txDataBits(1:min(length(rxDataBits_eve), length(txDataBits))));
            errors_eve = errors_eve + num_errors_eve;
            total_bits_eve = total_bits_eve + min(length(rxDataBits_eve), length(txDataBits));
        end

        BER_no_comp(fd_idx, snr_idx) = errors_no_comp / total_bits_no_comp;
        BER_with_comp(fd_idx, snr_idx) = errors_with_comp / total_bits_with_comp;
        BER_with_kalman(fd_idx, snr_idx) = errors_with_kalman / total_bits_with_kalman;
        BER_eve(fd_idx, snr_idx) = errors_eve / total_bits_eve;
        
        actual_BER_no_comp(fd_idx, snr_idx) = errors_no_comp / total_bits_no_comp;
        actual_BER_with_comp(fd_idx, snr_idx) = errors_with_comp / total_bits_with_comp;
        actual_BER_with_kalman(fd_idx, snr_idx) = errors_with_kalman / total_bits_with_kalman;
        actual_BER_eve(fd_idx, snr_idx) = errors_eve / total_bits_eve;
        
        BER_no_comp(fd_idx, snr_idx) = max(errors_no_comp / total_bits_no_comp, 1e-7);
        BER_with_comp(fd_idx, snr_idx) = max(errors_with_comp / total_bits_with_comp, 1e-7);
        BER_with_kalman(fd_idx, snr_idx) = max(errors_with_kalman / total_bits_with_kalman, 1e-7);
        BER_eve(fd_idx, snr_idx) = max(errors_eve / total_bits_eve, 1e-7);

        % Compute average phase MSE
        if phase_samples > 0
            phase_mse_basic(fd_idx, snr_idx) = phase_errors_basic / phase_samples;
            phase_mse_kalman(fd_idx, snr_idx) = phase_errors_kalman / phase_samples;
        end

        fprintf('No Comp BER = %.6e, Basic Comp BER = %.6e, Kalman Comp BER = %.6e, Eve BER = %.6e, Phase MSE (Basic/Kalman) = %.6e/%.6e\n', ...
                BER_no_comp(fd_idx, snr_idx), BER_with_comp(fd_idx, snr_idx), ...
                BER_with_kalman(fd_idx, snr_idx), BER_eve(fd_idx, snr_idx), ...
                phase_mse_basic(fd_idx, snr_idx), phase_mse_kalman(fd_idx, snr_idx));
    end
    fprintf('\n');
end

%% Plot Results
colors = {'k', 'r', 'g', 'b', 'm'};
markers = {'o', 's', 'd', '^', 'v'};

% DOA Estimation Error
figure;
for fd_idx = 1:length(doppler_fd_levels)
    plot(SNR_dB, doa_errors(fd_idx, :), [colors{fd_idx}, '-', markers{fd_idx}], 'LineWidth', 2, ...
        'DisplayName', ['fd = ', num2str(doppler_fd_levels(fd_idx))]);
    hold on;
end
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Mean Angular Error (degrees)', 'FontSize', 12);
title('DOA Estimation Error vs SNR for Different Doppler Frequencies', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10);

% BER Performance
for fd_idx = 1:length(doppler_fd_levels)
    figure;
    BER_no_comp_plot = BER_no_comp(fd_idx,:);
    BER_with_comp_plot = BER_with_comp(fd_idx,:);
    BER_with_kalman_plot = BER_with_kalman(fd_idx,:);
    BER_eve_plot = BER_eve(fd_idx,:);
    
    semilogy(SNR_dB, BER_no_comp_plot, 'r-^', 'LineWidth', 2, 'DisplayName', 'Without Compensation');
    hold on;
    semilogy(SNR_dB, BER_with_comp_plot, 'g-o', 'LineWidth', 2, 'DisplayName', 'Basic Compensation');
    semilogy(SNR_dB, BER_with_kalman_plot, 'b-s', 'LineWidth', 2, 'DisplayName', 'Kalman Compensation');
    semilogy(SNR_dB, BER_eve_plot, 'k-x', 'LineWidth', 2, 'DisplayName', 'Eavesdropper');
    for i = 1:length(SNR_dB)
        if actual_BER_with_comp(fd_idx, i) == 0
            plot(SNR_dB(i), 1e-7, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
        end
        if actual_BER_with_kalman(fd_idx, i) == 0
            plot(SNR_dB(i), 1e-7, 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
        end
    end
    grid on;
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('Bit Error Rate (BER)', 'FontSize', 12);
    title(['MIMO BER Performance with AN - Doppler fd = ', num2str(doppler_fd_levels(fd_idx))], 'FontSize', 14);
    legend('Location', 'southwest', 'FontSize', 10);
    axis([min(SNR_dB) max(SNR_dB) 1e-7 1]);
end

% Secrecy Gap
figure;
for fd_idx = 2:length(doppler_fd_levels)
    secrecy_gap = BER_eve(fd_idx,:) - BER_with_kalman(fd_idx,:);
    plot(SNR_dB, secrecy_gap, [colors{fd_idx}, '-', markers{fd_idx}], 'LineWidth', 2, ...
        'DisplayName', ['fd = ', num2str(doppler_fd_levels(fd_idx))]);
    hold on;
end
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Secrecy Gap (Eve BER - Legitimate BER)', 'FontSize', 12);
title('Security Performance with Artificial Noise', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);

% Combined Kalman Plot
figure;
semilogy(SNR_dB, BER_no_comp(1,:), 'k-o', 'LineWidth', 2, 'DisplayName', 'No Doppler');
hold on;
for fd_idx = 2:length(doppler_fd_levels)
    semilogy(SNR_dB, BER_with_kalman(fd_idx,:), [colors{fd_idx}, '-', markers{fd_idx}], 'LineWidth', 2, ...
        'DisplayName', ['fd = ', num2str(doppler_fd_levels(fd_idx)), ' (Kalman)']);
end
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('MIMO BER Performance with Kalman-based Doppler Compensation', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 10);
axis([min(SNR_dB) max(SNR_dB) 1e-7 1]);

% Secrecy Capacity
figure;
secrecy_capacity = zeros(length(doppler_fd_levels), length(SNR_dB));
for fd_idx = 1:length(doppler_fd_levels)
    for snr_idx = 1:length(SNR_dB)
        if BER_with_kalman(fd_idx, snr_idx) < 1e-6
            C_legitimate = log2(modulation_order);
        else
            C_legitimate = log2(modulation_order) * (1 - BER_with_kalman(fd_idx, snr_idx));
        end
        C_eve = log2(modulation_order) * (1 - BER_eve(fd_idx, snr_idx));
        secrecy_capacity(fd_idx, snr_idx) = max(0, C_legitimate - C_eve);
    end
    plot(SNR_dB, secrecy_capacity(fd_idx, :), [colors{fd_idx}, '-', markers{fd_idx}], 'LineWidth', 2, ...
        'DisplayName', ['fd = ', num2str(doppler_fd_levels(fd_idx))]);
    hold on;
end
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Secrecy Capacity (bits/symbol)', 'FontSize', 12);
title('Secrecy Capacity vs SNR with Artificial Noise', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
ylim([0 log2(modulation_order)]);

% Phase MSE Plot
figure;
for fd_idx = 1:length(doppler_fd_levels)
    if doppler_fd_levels(fd_idx) > 0
        semilogy(SNR_dB, phase_mse_basic(fd_idx, :), [colors{fd_idx}, '--'], 'LineWidth', 2, ...
            'DisplayName', ['Basic, fd = ', num2str(doppler_fd_levels(fd_idx))]);
        hold on;
        semilogy(SNR_dB, phase_mse_kalman(fd_idx, :), [colors{fd_idx}, '-'], 'LineWidth', 2, ...
            'DisplayName', ['Kalman, fd = ', num2str(doppler_fd_levels(fd_idx))]);
    end
end
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Phase MSE (radians^2)', 'FontSize', 12);
title('Phase Tracking Accuracy vs SNR', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10);

%% Improved DOA Estimation Function (MUSIC Algorithm)
function [est_angles, P_MUSIC] = improved_doa_estimation(rx_signal, antenna_spacing, n_angles, angle_range_rad, Ns)
    % Compute covariance matrix
    R = (rx_signal * rx_signal') / size(rx_signal, 2);
    
    % Eigenvalue decomposition
    [V, D] = eig(R);
    [~, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    
    % Noise subspace
    En = V(:, Ns+1:end);
    
    % MUSIC spectrum
    P_MUSIC = zeros(1, n_angles);
    for i = 1:n_angles
        theta = angle_range_rad(i);
        a = exp(-1j * 2 * pi * antenna_spacing * (0:size(rx_signal, 1)-1)' * cos(theta));
        P_MUSIC(i) = 1 / abs(a' * (En * En') * a);
    end
    
    % Normalize spectrum
    P_MUSIC = P_MUSIC / max(P_MUSIC);
    
    % Find peaks
    [~, locs] = findpeaks(P_MUSIC, 'SortStr', 'descend', 'NPeaks', 2);
    if length(locs) < 2
        est_angles = [rad2deg(angle_range_rad(locs(1))), 0];
    else
        est_angles = sort(rad2deg(angle_range_rad(locs)));
    end
end
