% MATLAB Simulation for Secure Wireless Communication System with Directional Modulation
% Thesis: Enhancing Physical Layer Security through Advanced Techniques for DOA Estimation
% Under Doppler Shift in Directional Modulation Systems
% Authors: Venkata Bhamidipati & Dheeraj Audam

clear; clc; close all;

%% === System Parameters ===
SNR_dB = 0:10:30;             % Reduced SNR points (1 x 4 array)
SNR = 10.^(SNR_dB/10);        % Linear SNR
fc = 2e9;                     % Carrier frequency (2 GHz)
v = 30;                       % Velocity (m/s)
c = 3e8;                      % Speed of light
lambda = c / fc;              % Wavelength
d = lambda / 2;               % Antenna spacing
Fs = 1e6;                     % Sampling frequency (1 MHz)
Ts = 1 / Fs;                  % Sampling period
nAntennas = 8;                % Number of antennas
true_doa_bob = 30;            % True DOA for Bob (degrees)
true_doa_eve = -20;           % True DOA for Eve (degrees)
f_d_bob = v * fc / c;         % Doppler shift for Bob (Hz)
f_d_eve = v * fc / c * 0.5;   % Doppler shift for Eve (Hz)
angles = -90:2:90;            % Coarser angle grid for MUSIC
K = 648;                      % Symbols to process (subset of full codeword)
P_s = 1;                      % Total transmit power
beta1 = sqrt(0.9);            % Power allocation for message
beta2 = sqrt(0.1);            % Power allocation for AN

%% LDPC Setup (DVB-S.2, Rate 1/2, Full Matrix)
H = dvbs2ldpc(1/2, 'sparse'); % 32400 x 64800 matrix
N_ldpc = 64800;               % Full codeword length (bits)
K_ldpc = 32400;               % Full message length (bits)
ldpcEncCfg = ldpcEncoderConfig(H);
ldpcDecCfg = ldpcDecoderConfig(H);
maxIter = 50;                 % Reduced iterations for decoding

%% Simulation Settings
nTrials = 5;                  % Reduced trials
numSnapshots = 100;           % Reduced snapshots
M = nAntennas;
ber_bob = zeros(size(SNR));
ber_eve = zeros(size(SNR));
rmse_theta_bob = zeros(size(SNR));
crlb_theta_bob = zeros(size(SNR));

%% Steering Vector Function
steerVec = @(theta) exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta) / lambda) / sqrt(M);

%% QPSK Setup
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

%% Precompute for Speed
t = (0:K-1)' * Ts;            % K x 1 (648 symbols)

for idx = 1:length(SNR_dB)
    fprintf('Processing SNR = %d dB (%d of %d)\n', SNR_dB(idx), idx, length(SNR_dB));
    snr = SNR(idx);
    noise_var = P_s / snr;
    ber_bob_sum = 0;
    ber_eve_sum = 0;
    doa_errors_bob = zeros(1, nTrials);

    for trial = 1:nTrials
        %% Transmitter
        bits = randi([0 1], K_ldpc, 1); % 32400 x 1 (full message)
        encoded = ldpcEncode(bits, ldpcEncCfg); % 64800 x 1 (full codeword)
        modData_full = qpskMod(encoded); % 32400 x 1 (N_ldpc / 2 symbols)
        modData = modData_full(1:K); % 648 x 1 (subset for transmission)

        %% DOA & Doppler Estimation
        s_bob = modData(1:min(numSnapshots, K)); % 100 x 1
        s_eve = qpskMod(randi([0 1], numSnapshots * 2, 1)); % 100 x 1
        Y = zeros(M, numSnapshots);
        h_bob = steerVec(true_doa_bob); % 8 x 1
        h_eve = steerVec(true_doa_eve); % 8 x 1
        alpha_b = sqrt(0.5) * (randn(1, numSnapshots) + 1j * randn(1, numSnapshots));
        alpha_e = sqrt(0.5) * (randn(1, numSnapshots) + 1j * randn(1, numSnapshots));
        doppler_b = exp(1j * 2 * pi * f_d_bob * (0:numSnapshots-1) * Ts); % 1 x 100
        doppler_e = exp(1j * 2 * pi * f_d_eve * (0:numSnapshots-1) * Ts); % 1 x 100
        Y = h_bob * (alpha_b .* s_bob.' .* doppler_b) + ...
            h_eve * (alpha_e .* s_eve.' .* doppler_e) + ...
            sqrt(noise_var / 2) * (randn(M, numSnapshots) + 1j * randn(M, numSnapshots)); % 8 x 100

        % MUSIC DOA (Vectorized)
        Ryy = (Y * Y') / numSnapshots; % 8 x 8
        [V, D] = eig(Ryy);
        [~, idxs] = sort(diag(D), 'descend');
        Vn = V(:, idxs(3:end)); % 8 x 6
        A = exp(-1j * 2 * pi * d * (0:M-1)' * sind(angles) / lambda) / sqrt(M); % 8 x 91
        Pmusic = 1 ./ sum(abs(A' * Vn).^2, 2)'; % 1 x 91
        [~, peaks] = findpeaks(Pmusic, 'SortStr', 'descend', 'NPeaks', 2);
        est_doas = angles(peaks);
        dist_to_bob = abs(est_doas - true_doa_bob);
        [~, idx_b] = min(dist_to_bob);
        est_doa_bob = est_doas(idx_b);
        est_doa_eve = est_doas(3 - idx_b);
        doa_errors_bob(trial) = est_doa_bob - true_doa_bob;

        % CRLB
        Lambda = sort(diag(D), 'descend');
        sigma_hat2 = mean(Lambda(3:end));
        gamma_b = (Lambda(1) - sigma_hat2) / (M * sigma_hat2);
        d_n_sq_sum = sum((((1:M) - (M + 1)/2) * d).^2);
        crlb_theta_bob(idx) = (lambda^2 / (8 * pi^2 * numSnapshots * gamma_b * cosd(est_doa_bob)^2 * d_n_sq_sum));

        % Doppler Estimation (Vectorized)
        h_est_bob = steerVec(est_doa_bob); % 8 x 1
        s_est_bob = h_est_bob' * Y; % 1 x 100
        z_k_bob = s_est_bob ./ s_bob'; % 1 x 100
        f_range = linspace(-Fs/2, Fs/2, 101); % 101 x 1
        A_f = exp(1j * 2 * pi * f_range' * (0:numSnapshots-1) * Ts); % 101 x 100
        Pmusic_f = 1 ./ sum(abs(A_f * z_k_bob(:)).^2, 2)'; % 1 x 101
        [~, f_peak] = max(Pmusic_f);
        f_d_bob_est = f_range(f_peak);

        %% Beamforming (Simplified)
        delta_theta_max_bob = sqrt(crlb_theta_bob(idx));
        delta_theta_max_eve = delta_theta_max_bob;
        u_b = steerVec(est_doa_bob); % 8 x 1
        u_e = steerVec(est_doa_eve); % 8 x 1
        R_e = eye(M) - u_e * u_e'; % 8 x 8
        R_b = eye(M) - u_b * u_b'; % 8 x 8
        vb = R_e * u_b / norm(R_e * u_b); % 8 x 1
        we = R_b * u_e / norm(R_b * u_e); % 8 x 1

        %% Transmit Signal
        s_data = modData; % 648 x 1
        an_noise = (randn(K, 1) + 1j * randn(K, 1)) / sqrt(2); % 648 x 1
        txSig = beta1 * sqrt(P_s) * vb * s_data.' + beta2 * sqrt(P_s) * we * an_noise.'; % 8 x 648

        %% Channel Effect
        channel_bob = h_bob' * txSig; % 1 x 648
        doppler_bob = exp(1j * 2 * pi * f_d_bob * t); % 648 x 1
        rxSig_bob = alpha_b(1) * (channel_bob.' .* doppler_bob) + ...
                    sqrt(noise_var / 2) * (randn(K, 1) + 1j * randn(K, 1)); % 648 x 1
        channel_eve = h_eve' * txSig; % 1 x 648
        doppler_eve = exp(1j * 2 * pi * f_d_eve * t); % 648 x 1
        rxSig_eve = alpha_e(1) * (channel_eve.' .* doppler_eve) + ...
                    sqrt(noise_var / 2) * (randn(K, 1) + 1j * randn(K, 1)); % 648 x 1

        %% Receiver Processing
        doppler_comp_bob = exp(-1j * 2 * pi * f_d_bob_est * t); % 648 x 1
        rxComp_bob = rxSig_bob .* doppler_comp_bob; % 648 x 1
        rxBits_bob = qpskDemod(rxComp_bob); % 1296 x 1 (2 bits per symbol)
        % Pad with zeros to match N_ldpc for decoding
        rxBits_bob_full = [rxBits_bob; zeros(N_ldpc - length(rxBits_bob), 1)]; % 64800 x 1
        decoded_bob = ldpcDecode(rxBits_bob_full, ldpcDecCfg, maxIter); % 32400 x 1
        decoded_bob_subset = decoded_bob(1:K_ldpc); % 32400 x 1 (full message)

        doppler_comp_eve = exp(-1j * 2 * pi * f_d_eve * t); % 648 x 1
        rxComp_eve = rxSig_eve .* doppler_comp_eve; % 648 x 1
        rxBits_eve = qpskDemod(rxComp_eve); % 1296 x 1
        rxBits_eve_full = [rxBits_eve; zeros(N_ldpc - length(rxBits_eve), 1)]; % 64800 x 1
        decoded_eve = ldpcDecode(rxBits_eve_full, ldpcDecCfg, maxIter); % 32400 x 1
        decoded_eve_subset = decoded_eve(1:K_ldpc); % 32400 x 1

        % BER Calculation (on subset of bits)
        ber_bob_sum = ber_bob_sum + biterr(decoded_bob_subset(1:K/2), bits(1:K/2)) / (K/2);
        ber_eve_sum = ber_eve_sum + biterr(decoded_eve_subset(1:K/2), bits(1:K/2)) / (K/2);
    end

    ber_bob(idx) = ber_bob_sum / nTrials;
    ber_eve(idx) = ber_eve_sum / nTrials;
    rmse_theta_bob(idx) = sqrt(mean(doa_errors_bob.^2));
end

% Check SNR_dB after loop
if ~exist('SNR_dB', 'var') || isempty(SNR_dB)
    error('SNR_dB lost after simulation loop. It was cleared unexpectedly.');
end
fprintf('Simulation complete. SNR_dB = %s\n', mat2str(SNR_dB));

%% === PLOTS ===
% Debug Check for SNR_dB
if ~exist('SNR_dB', 'var') || isempty(SNR_dB)
    error('SNR_dB is not defined or is empty before plotting. Check if it was cleared.');
end
if ~isvector(SNR_dB) || length(SNR_dB) ~= length(ber_bob)
    error('SNR_dB size mismatch. Expected length %d, got %s', length(ber_bob), mat2str(size(SNR_dB)));
end

figure('Name', 'AR-NSP Secure Transmission Performance', 'NumberTitle', 'off');
subplot(1, 3, 1);
semilogy(SNR_dB, ber_bob, 'bo-', 'LineWidth', 2); hold on;
semilogy(SNR_dB, ber_eve, 'rx--', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('Bit Error Rate');
title('BER vs SNR (AR-NSP)');
legend('Bob', 'Eve'); grid on;

subplot(1, 3, 2);
plot(SNR_dB, rmse_theta_bob, 'rs--', 'LineWidth', 2); hold on;
plot(SNR_dB, sqrt(crlb_theta_bob), 'k-.', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('DOA RMSE (deg)');
title('DOA Estimation RMSE vs SNR (Bob)');
legend('RMSE', 'CRLB'); grid on;

subplot(1, 3, 3);
plot(SNR_dB, sqrt(crlb_theta_bob), 'g^-', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('CRLB (deg)');
title('CRLB vs SNR (Bob)'); grid on;

%% Function Definitions (Simplified)
function R = compute_R(~, ~, M, ~, ~)
    R = eye(M); % Simplified for speed
end

function u = compute_u(theta_hat, ~, M, d, lambda)
    u = exp(-1j * 2 * pi * d * (0:M-1)' * sind(theta_hat) / lambda) / sqrt(M);
end
