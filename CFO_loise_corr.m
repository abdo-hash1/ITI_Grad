% --- QPSK Transceiver with Coarse CFO ---

clear; clc; 

% --- 1. System Parameters ---
Rs = 1.92e6 ;           % Symbol rate
Ts = 1/Rs;
alpha = 0.5;            % Roll-off factor for RRC filter
N_span = 6;             % Filter span (symbols)
OSR = 4;                % Oversampling factor
N_payload_bits = 2240;  % Fixed payload size in bits

SNR_dB_range = 0:1:60;  

% ---- CFO PARAMETER (using PPM) ----
% noise = randn([1,4608]);
% Choose a nominal RF carrier frequency for your system (in Hz)
fc_hz = 2.4e9; 

% Define your desired ppm value
ppm_offset = 25;

% Calculate the absolute frequency offset in Hz
max_offset_hz = (ppm_offset * 1e-6) * fc_hz;
step_hz = 1;    

% Create the vector of frequency offsets
cfo_range_hz = -max_offset_hz : step_hz : max_offset_hz;




BER = zeros(size(SNR_dB_range));
BER_theoretical = zeros(size(SNR_dB_range));



% --- Preamble, Guard, Filter Design --- (Same as your working version)
% Define Barker sequence using binary 0s and 1s
barker = [1 1 1 1 1 0 0 1 1 0 1 0 1];

% Repeat each bit twice
repeated_bits = repelem(barker, 2);

% Convert to bipolar [-1, 1]
bipolar_seq = (2 * repeated_bits - 1); 

% I/Q Mapping
unnormalized_symbols = bipolar_seq + 1j * bipolar_seq;

% Power Normalization
preamble_symbols_complex = unnormalized_symbols / sqrt(2);

num_preamble_symbols = length(preamble_symbols_complex); 
fprintf('Preamble generated: %d symbols.\n', num_preamble_symbols);

% Separate I and Q components if needed by your script structure
I_preamble_sym = real(preamble_symbols_complex);
Q_preamble_sym = imag(preamble_symbols_complex);

[h_filter, ~] = designPulseShapingFilter(Ts, alpha, N_span, OSR, 'RRC');
total_filter_delay_samples = length(h_filter) - 1;





%% Main Simulation Loop
for j_idx = 1:length(SNR_dB_range)
    random_cfo_value = cfo_range_hz(1000);
    current_SNR = SNR_dB_range(j_idx); 
    fprintf('\n\nSimulating for SNR = %.1f dB and CFO =%.1f ...\n', current_SNR,random_cfo_value);
    % ==== Transmit Path ==== 
    % Transmitted Bits
    payload_data_tx_bits = randi([0 1], 1, N_payload_bits);
    BitsPerSymbol = 2;
    % I/Q Modulation
    [I_payload, Q_payload] = QPSK(payload_data_tx_bits);
    % Adding Preamble to Payload 
    payload_symbols_qpsk_complex = I_payload + 1j*Q_payload;
    num_payload_symbols = length(payload_symbols_qpsk_complex);
    signal_power=var(payload_symbols_qpsk_complex);
     % fprintf('\n\nSignal Power after Mod =%.1f ...\n', signal_power);

    total_symbols_complex_tx = [preamble_symbols_complex,payload_symbols_qpsk_complex];
    signal_power=var(total_symbols_complex_tx);
    % fprintf('\n\nSignal Power after adding preamble =%.1f ...\n', signal_power);

    % Upsampling
    upsampled_total_symbols_complex = upsample(total_symbols_complex_tx, OSR);
    % upsampled_total_symbols_complex = upsampled_total_symbols_complex .*sqrt(OSR);
    signal_power=var(upsampled_total_symbols_complex);
    % fprintf('\n\nSignal Power after upsampling =%.10f ...\n', signal_power);

    % Pulse Shaping
    tx_signal_shaped = conv(upsampled_total_symbols_complex, h_filter);
    signal_power=var(tx_signal_shaped);
    % fprintf('\n\nSignal Power after pulse Shapping =%.10f ...\n', signal_power);
    signal_power=var(tx_signal_shaped*OSR);
    % fprintf('\n\nSignal Power after  Filtering scale =%.10f ...\n', signal_power);
    tx_signal_shaped = tx_signal_shaped.*OSR;
    % ==== Channel Model ==== 
    % Adding CFO
    num_tx_samples = length(tx_signal_shaped);
    t_osr_tx = (0:num_tx_samples-1) * (Ts/OSR);
    tx_signal_with_cfo = tx_signal_shaped .* exp(1j * 2 * pi * random_cfo_value * t_osr_tx);
    % Adding AWGN Noise
    
    %snr_for_awgn_func_dB = current_SNR + 10*log10(BitsPerSymbol / OSR);
    noisy_signal = adding_awgn_noise(tx_signal_with_cfo,current_SNR) ;
    signal_power=var(noisy_signal);
    % fprintf('\n\nSignal Power after Noise Adding =%.10f ...\n', signal_power);
    % ==== RECEIVER ====
    rx_signal_filtered_complex = conv(noisy_signal, h_filter);
    signal_power=var(rx_signal_filtered_complex);
    % fprintf('\n\nSignal Power after Matched Filtering =%.10f ...\n', signal_power);
  
    
    % rx_signal_filtered_complex = rx_signal_filtered_complex ./ sqrt(OSR) ;
    % Downsampling & Truncation
    total_delay = length(h_filter) - 1;
    rx_downsampled_long = rx_signal_filtered_complex(total_delay + 1 : OSR : end);
    num_tx_symbols = length(total_symbols_complex_tx);
    rx_downsampled = rx_downsampled_long(1:num_tx_symbols);
    signal_power=var(rx_downsampled);
    % fprintf('\n\nSignal Power after Down Sampling =%.1f ...\n', signal_power);

    % Coarse CFO Estimation 
    N = 26; % Number of Preamble
    M = N/2; % Number of Autocorrlation taps
    rx_preamble = rx_downsampled(1:N);
    % Remove modulation from Preamble Symbols
    x = rx_preamble .* conj(preamble_symbols_complex);
    
    R = zeros(M, 1);
    for k = 1:M
        temp_sum = 0;
        for i = k+1:N
         temp_sum = temp_sum + (x(i) * conj(x(i-k)));
        end
        R(k) = (1/(N-k)) * temp_sum;
    end
    S = sum(R);
    delta_f_hat = (1 / (pi * Ts * (M + 1))) * angle(S);
    fprintf('  Actual CFO: %.2f Hz | Manual Estimate: %.2f Hz | Difference: %.2f Hz  \n', random_cfo_value, delta_f_hat, abs(random_cfo_value-delta_f_hat));
    % fprintf('  Angle (S): %.10f Hz | CORDIC angle: %.10f Hz | Difference: %.10f Hz  \n', angle(S), deg2rad(CORDIC_Vector(S)), abs(angle(S)-deg2rad(CORDIC_Vector(S))));
    
    % Coarse CFO Correction (applied to the whole packet)
    t_full_packet = (0:length(rx_downsampled)-1).' * Ts;
    cfo_corrector_full = exp(-1j * 2 * pi * delta_f_hat * t_full_packet);
    rx_packet_coarse_corrected = rx_downsampled(:) .* cfo_corrector_full(:);
    
    
    %  Extract the Payload for Demodulation
    payload_symbols_for_demod = rx_packet_coarse_corrected(N+1:end);
    
    % --- Demodulation ---
    rx_payload_bits = QPSK_dem(real(payload_symbols_for_demod), imag(payload_symbols_for_demod));

    % ==== BER Calculation on PAYLOAD ====
    num_errors = sum(rx_payload_bits ~= payload_data_tx_bits);
    BER_simulated(j_idx) = num_errors / N_payload_bits;
    fprintf('  BER: %.10e\n', BER_simulated(j_idx));

    % ==== Theoretical BER ====
    ebno_linear = 10^(current_SNR/10);
    BER_theoretical(j_idx) = qfunc(sqrt(2 * ebno_linear));
end % End of SNR loop

%% Plotting
% BER Plot
figure;
semilogy(SNR_dB_range, BER_simulated, 'bo-', 'LineWidth', 1.5);
hold on; semilogy(SNR_dB_range, BER_theoretical, 'r--', 'LineWidth', 1.5);
hold off; grid on; xlabel('SNR (dB)'); ylabel('BER');
title(sprintf('QPSK BER with CFO Correction'));
legend('Simulated BER', 'Theoretical BER', 'Location','southwest');
ylim_min_ber = max(1e-6, min(BER_simulated(BER_simulated>0))/2); if isempty(ylim_min_ber), ylim_min_ber = 1e-6; end; ylim([ylim_min_ber 1]);
