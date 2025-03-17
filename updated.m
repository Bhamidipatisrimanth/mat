clc; clear; close all;

%% PARAMETERS
numBits = 1024;        % Number of bits per message
M = 4;                 % QPSK
SNRdB = 0:5:30;        % SNR range
numAntennasTx = 8;     % Tx antennas
numAntennasRx = 8;     % Rx antennas (used for MUSIC)
numPaths = 3;          % Multipath components
velocity = 30;         % m/s
fc = 2e9;              % Carrier frequency
c = 3e8;               % Speed of light
alpha_AN = 0.3;        % Artificial Noise scaling
theta_target = 30;     % DOA target angle in degrees

f_D = velocity / c * fc;   % Doppler shift (Hz)

%% LDPC Parameters (Custom)
k = 8; n = 16;          % Simple LDPC for demonstration
H = ldpcParityCheckMatrix(k, n);

%% Initialize BER Arrays
ber_legit_results = zeros(length(SNRdB), 1);
ber_eve_results = zeros(length(SNRdB), 1);

%% MAIN LOOP OVER SNR VALUES
for idx = 1:length(SNRdB)
    SNR = SNRdB(idx);

    %% 1. Data Generation & LDPC Encoding
    data = randi([0 1], k, 1);
    codedData = ldpcEncodeCustom(data, H);

    %% 2. QPSK Modulation
    bitsPerSymbol = log2(M);
    symbolsMatrix = reshape([codedData; zeros(mod(-numel(codedData), bitsPerSymbol), 1)], bitsPerSymbol, []).';
    modulatedSymbols = pskmod(bi2de(symbolsMatrix, 'left-msb'), M, pi/M);

    %% 3. MIMO Precoding (MRT)
    W = randn(numAntennasTx, 1) + 1j * randn(numAntennasTx, 1);
    W = W / norm(W);  % Normalize
    bfSignal = W * modulatedSymbols.'; % [Tx antennas x NumSymbols]

    %% 4. Directional Modulation (DM)
    steeringVec = exp(1j * pi * (0:numAntennasTx-1).' * sind(theta_target));
    dmSignal = bfSignal .* steeringVec; % Element-wise multiplication for directional modulation

    %% 5. Artificial Noise (AN)
    nullSpace = null(W.');
    AN = alpha_AN * nullSpace * (randn(size(nullSpace,2), length(modulatedSymbols)) + 1j * randn(size(nullSpace,2), length(modulatedSymbols)));

    %% 6. Transmit Signal (DM + AN)
    txSignal = dmSignal + AN; % Ensure correct matrix dimensions

    %% 7. Channel Model (Multipath + Doppler)
    h = (randn(numPaths, 1) + 1j * randn(numPaths, 1)) / sqrt(2);
    t = (0:length(modulatedSymbols)-1).';
    rxSignal = zeros(size(modulatedSymbols));
    for p = 1:numPaths
        dopplerShift = exp(1j * 2 * pi * f_D * t);
        delayedSignal = [zeros(p, 1); txSignal(1, 1:end-p).'];
        rxSignal = rxSignal + h(p) * delayedSignal .* dopplerShift;
    end

    %% 8. AWGN
    noisePower = norm(rxSignal)^2 / (length(rxSignal) * 10^(SNR/10));
    noise = sqrt(noisePower/2) * (randn(size(rxSignal)) + 1j * randn(size(rxSignal)));
    rx_noisy = rxSignal + noise;

    %% 9. MUSIC DOA Estimation + Kalman Filter
    [angles, spectrum] = musicDOA(rx_noisy, numAntennasRx, numPaths);
    smoothedAngles = kalmanFilterDOA(angles);

    %% 10. Doppler Estimation + Kalman Filter
    dopplerEst = kalmanDoppler(rx_noisy, length(rx_noisy));
    rxComp = rx_noisy .* exp(-1j * 2 * pi * dopplerEst * t);

    %% 11. Zero Forcing Equalization (Fixed Dimension Error)
    eqSignal_legit = pinv(W) * rxComp; % Ensure matrix dimensions match

    %% 12. Demodulation & LDPC Decoding
    demodSymbols = pskdemod(eqSignal_legit, M, pi/M);
    demodBits = reshape(de2bi(demodSymbols, bitsPerSymbol, 'left-msb').', [], 1);
    decodedData_legit = ldpcDecodeCustom(demodBits(1:k), H);

    %% 13. BER Calculation
    ber_legit_results(idx) = sum(decodedData_legit ~= data) / length(data);
    fprintf('SNR: %d dB | Legit BER: %.4e\n', SNR, ber_legit_results(idx));
end

%% 14. BER vs. SNR Plot (Fix visibility issue)
figure;
semilogy(SNRdB, ber_legit_results, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Legitimate Receiver');
xlabel('SNR (dB)'); ylabel('BER'); title('BER vs. SNR'); legend show; grid on;
ylim([1e-5, 1]); % Ensure proper scaling
drawnow; % Force figure rendering

%% 15. Constellation Diagrams (Ensure correct data is plotted)
figure;
scatterplot(eqSignal_legit(:));
title('Constellation: Legitimate Receiver');
drawnow;
pause(0.1);

figure;
scatterplot(rx_noisy(:));
title('Constellation: Eavesdropper');
drawnow;
pause(0.1);

%% 16. MUSIC DOA Spectrum Plot (Fix range and visibility)
figure;
plot(angles, 10*log10(abs(spectrum)), 'LineWidth', 2);
xlabel('Angle (Degrees)');
ylabel('Pseudo Spectrum (dB)');
title('MUSIC DOA Spectrum');
grid on;
drawnow;
