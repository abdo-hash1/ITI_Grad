% --- QPSK Transceiver with Coarse CFO---

clear; clc; 

% --- 1. System Parameters ---
Rs = 1.92e6 ;           % Symbol rate
Ts = 1/Rs;
alpha = 0.5;            % Roll-off factor for RRC filter
N_span = 6;             % Filter span (symbols)
OSR = 4;                % Oversampling factor
N_payload_bits = 2240;  % Fixed payload size in bits

SNR_dB_range = 0:1:60;  % SNR range for BER curve (e.g., Eb/No)

best_known_sample_offset = 0; 

BER = zeros(size(SNR_dB_range));
BER_theoretical = zeros(size(SNR_dB_range));



% --- Preamble, Filter Design --- 
% 1. Define Barker sequence using binary 0s and 1s
barker = [1 1 1 1 1 0 0 1 1 0 1 0 1];

% 2. Repeat each bit twice
repeated_bits = repelem(barker, 2);

% 3. Convert binary [0, 1] to bipolar [-1, 1]
bipolar_seq = (2 * repeated_bits - 1); % Keep as row vector for now

% 4. Apply the I=Q mapping rule
unnormalized_symbols = bipolar_seq + 1j * bipolar_seq;

% 5. Normalize the symbols and store in the final variable
preamble_symbols_complex = unnormalized_symbols / sqrt(2);

num_preamble_symbols = length(preamble_symbols_complex); % Should be 26
fprintf('Preamble generated: %d symbols.\n', num_preamble_symbols);

% Separate I and Q components if needed by your script structure
I_preamble_sym = real(preamble_symbols_complex);
Q_preamble_sym = imag(preamble_symbols_complex);

[h_filter, ~] = designPulseShapingFilter(Ts, alpha, N_span, OSR, 'RRC');
total_filter_delay_samples = length(h_filter) - 1;


%% Main Simulation Loop
for j_idx = 1:length(SNR_dB_range)
    current_EbNo_dB = SNR_dB_range(j_idx); 
    fprintf('\nSimulating for Eb/No = %.1f dB...\n', current_EbNo_dB);
    % ==== Transmit Path ==== 
    payload_data_tx_bits = randi([0 1], 1, N_payload_bits);
    BitsPerSymbol = 2;
    refConstellation_qpsk_data = (1/sqrt(2)) * [1+1i, -1+1i, -1-1i, 1-1i];
    [I_payload, Q_payload] = QPSK(payload_data_tx_bits);
    
    payload_symbols_qpsk_complex = I_payload + 1j*Q_payload;
    num_payload_symbols = length(payload_symbols_qpsk_complex);
    total_symbols_complex_tx = [preamble_symbols_complex,payload_symbols_qpsk_complex];
    upsampled_total_symbols_complex = upsample(total_symbols_complex_tx, OSR);
    tx_signal_shaped = conv(upsampled_total_symbols_complex, h_filter);

    % ==== Channel Model ==== 

     % Adding AWGN Noise
    snr_for_awgn_func_dB = current_EbNo_dB + 10*log10(BitsPerSymbol / OSR);
    noisy_signal = awgn(tx_signal_shaped, snr_for_awgn_func_dB, 'measured');

    % ==== RECEIVER ====
    % rx_signal_filtered_complex = filter(h_filter, 1, noisy_signal);
    rx_signal_filtered_complex = conv(noisy_signal, h_filter);

    % ==== Down Sampling ====
    % --- Corrected Down Sampling and ROBUST TRUNCATION ---

    % 1. Calculate the total delay from BOTH Tx and Rx filters
    total_delay = length(h_filter) - 1;
    
    % 2. Perform the downsampling starting at the correct delay
    %    This vector might be slightly too long.
    rx_downsampled_long = rx_signal_filtered_complex(total_delay + 1 : OSR : end);
    
    % 3. **THIS IS THE KEY FIX:** Truncate the vector to the exact number of
    %    symbols that were transmitted.
    num_tx_symbols = length(total_symbols_complex_tx);
    rx_downsampled = rx_downsampled_long(1:num_tx_symbols);



    % rx_downsampled = rx_signal_filtered_complex(1:OSR:end);
    
    rx_payload = rx_downsampled(27:end) ;

   




    % ==== Demodulate PAYLOAD (using the final PLL-corrected symbols) ====
    rx_payload_bits = QPSK_dem(real(rx_payload), imag(rx_payload));
   
    % ==== BER Calculation on PAYLOAD ====
    num_errors = sum(rx_payload_bits ~= payload_data_tx_bits);
    BER_simulated(j_idx) = num_errors / N_payload_bits;
    fprintf('  BER: %.2e\n', BER_simulated(j_idx));

    % ==== Theoretical BER ====
    ebno_linear = 10^(current_EbNo_dB/10);
    BER_theoretical(j_idx) = qfunc(sqrt(2 * ebno_linear));
end % End of SNR loop

%% Plotting
% BER Plot
figure;
semilogy(SNR_dB_range, BER_simulated, 'bo-', 'LineWidth', 1.5);
hold on; semilogy(SNR_dB_range, BER_theoretical, 'r--', 'LineWidth', 1.5);
hold off; grid on; xlabel('Eb/No (dB)'); ylabel('BER');
title(sprintf('QPSK BER with AWGN only Channel'));
legend('Simulated BER', 'Theoretical BER', 'Location','southwest');
ylim_min_ber = max(1e-6, min(BER_simulated(BER_simulated>0))/2); if isempty(ylim_min_ber), ylim_min_ber = 1e-6; end; ylim([ylim_min_ber 1]);

