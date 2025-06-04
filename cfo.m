% --- Main QPSK Transceiver BER Simulation with CFO Estimation & Correction ---
% --- Includes Plot for Estimated CFO vs. SNR ---

clear; clc; close all;

% --- 1. System Parameters ---
Rs = 1.92e6 ;           % Symbol rate
Ts = 1/Rs;
alpha = 0.5;            % Roll-off factor for RRC filter
N_span = 6;             % Filter span (symbols) << Used for RRC filter & guard interval
OSR = 4;                % Oversampling factor
N_payload_bits = 2240;  % Fixed payload size in bits

SNR_dB_range = 0:2:60;  % SNR range for BER curve (e.g., Eb/No)

% --- USER-DEFINED SETTINGS based on previous debugging/choices ---
best_known_sample_offset = 0; % From your successful debugging

% Define the actual CFO to be introduced in the channel FOR THIS RUN
delta_f_actual_setting = 0.001 * Rs; % Example: 0.5% of Rs. Set to 0 to see estimation noise.
% delta_f_actual_setting = 0; 

enable_cfo_correction = true;    % Keep true to see corrected performance
% enable_cfo_correction = false; % To see effect of CFO without correction (est. CFO will be forced to 0)


fprintf('--- QPSK System BER Simulation ---\n');
fprintf('Symbol Rate: %.2f Msps, OSR: %d\n', Rs/1e6, OSR);
fprintf('Downsampling Offset: %d\n', best_known_sample_offset);
fprintf('Actual CFO in Channel: %.2f Hz (%.3f * Rs)\n', delta_f_actual_setting, delta_f_actual_setting/Rs);
if enable_cfo_correction
    fprintf('Receiver CFO Correction: ENABLED\n');
else
    fprintf('Receiver CFO Correction: DISABLED (Est. CFO forced to 0 for correction)\n');
end
fprintf('------------------------------------\n\n');

BER = zeros(size(SNR_dB_range));
evm_percent = zeros(size(SNR_dB_range));
estimated_cfo_values = zeros(size(SNR_dB_range)); % << NEW: To store estimated CFO

% --- Preamble Generation --- (Same as your working version)
barker_bits_map = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
K_preamble = 1/sqrt(2);
preamble_symbols_I = []; preamble_symbols_Q = [];
for bit_val = barker_bits_map
    symbol_I_comp = K_preamble * bit_val; symbol_Q_comp = K_preamble * bit_val;
    preamble_symbols_I = [preamble_symbols_I, symbol_I_comp, symbol_I_comp];
    preamble_symbols_Q = [preamble_symbols_Q, symbol_Q_comp, symbol_Q_comp];
end
preamble_symbols_complex = preamble_symbols_I + 1j*preamble_symbols_Q;
num_preamble_symbols = length(preamble_symbols_complex);

% --- Guard Symbols --- (Same as your working version)
num_guard_symbols = N_span; 
guard_symbols_complex = zeros(1, num_guard_symbols);

%% Design an RRC filter (Assuming designPulseShapingFilter.m is available)
if ~exist('designPulseShapingFilter','file'), error('designPulseShapingFilter.m not found.'); end
[h_filter, ~] = designPulseShapingFilter(Ts, alpha, N_span, OSR, 'RRC');
total_filter_delay_samples = length(h_filter) - 1;

%% Main Simulation Loop
for j_idx = 1:length(SNR_dB_range)
    current_EbNo_dB = SNR_dB_range(j_idx);
    fprintf('\nSimulating for Eb/No = %.1f dB...\n', current_EbNo_dB);

    % ==== Generate PAYLOAD Data and QPSK Symbols ====
    payload_data_bits = randi([0 1], 1, N_payload_bits);
    BitsPerSymbol_QPSK = 2;
    refConstellation_qpsk_data = (1/sqrt(2)) * [1+1i, -1+1i, -1-1i, 1-1i];
    if exist('QPSK_mod', 'file')
        [I_payload, Q_payload] = QPSK_mod(payload_data_bits);
    else % Fallback
        if j_idx==1, fprintf('FALLBACK: QPSK_mod.m not found. Using standard QPSK mapping for PAYLOAD.\n');end
        payload_symbols_qpsk_complex_fb = zeros(1, N_payload_bits/2);
        bit_pairs = reshape(payload_data_bits, 2, N_payload_bits/2)';
        for k_sym = 1:(N_payload_bits/2)
            if     isequal(bit_pairs(k_sym,:),[0 0]), payload_symbols_qpsk_complex_fb(k_sym) = refConstellation_qpsk_data(1);
            elseif isequal(bit_pairs(k_sym,:),[0 1]), payload_symbols_qpsk_complex_fb(k_sym) = refConstellation_qpsk_data(2);
            elseif isequal(bit_pairs(k_sym,:),[1 1]), payload_symbols_qpsk_complex_fb(k_sym) = refConstellation_qpsk_data(3);
            else, payload_symbols_qpsk_complex_fb(k_sym) = refConstellation_qpsk_data(4); end
        end
        I_payload = real(payload_symbols_qpsk_complex_fb); Q_payload = imag(payload_symbols_qpsk_complex_fb);
    end
    payload_symbols_qpsk_complex = I_payload + 1j*Q_payload;
    num_payload_symbols = length(payload_symbols_qpsk_complex);

    % ==== Combine, Upsample, Tx Filter ====
    total_symbols_complex_tx = [preamble_symbols_complex, guard_symbols_complex, payload_symbols_qpsk_complex];
    if ~exist('zero_stuffing','file')
        if j_idx==1, fprintf('Warning: zero_stuffing.m not found. Using MATLAB upsample().\n'); end
        upsampled_total_symbols_complex = upsample(total_symbols_complex_tx, OSR);
    else
        upsampled_total_symbols_complex = zero_stuffing(real(total_symbols_complex_tx), OSR) + ...
                                     1j*zero_stuffing(imag(total_symbols_complex_tx), OSR);
    end
    tx_signal_shaped = filter(h_filter, 1, upsampled_total_symbols_complex); % Using filter()

    % ==== CHANNEL MODEL ====
    num_tx_samples = length(tx_signal_shaped);
    t_osr_tx = (0:num_tx_samples-1) * (Ts/OSR);
    tx_signal_with_cfo = tx_signal_shaped .* exp(1j * 2 * pi * delta_f_actual_setting * t_osr_tx);

    EsNo_dB = current_EbNo_dB + 10*log10(BitsPerSymbol_QPSK);
    snr_for_awgn_func_dB = EsNo_dB - 10*log10(OSR); 
    if exist('adding_awgn_noise','file') && ~isempty(which('adding_awgn_noise'))
        noisy_signal = adding_awgn_noise(tx_signal_with_cfo, current_EbNo_dB);
    else
        if j_idx == 1, fprintf('Warning: adding_awgn_noise.m not found. Using MATLAB awgn() with calculated SNR.\n'); end
        noisy_signal = awgn(tx_signal_with_cfo, snr_for_awgn_func_dB, 'measured');
    end

    % ==== RECEIVER ====
    rx_signal_filtered_complex = filter(h_filter, 1, noisy_signal); % Using filter()

    % --- COARSE CFO ESTIMATION ---
    delta_f_estimated = 0; 
    if enable_cfo_correction
        start_idx_oversampled_preamble_in_rx_filt = 1 + total_filter_delay_samples;
        end_idx_oversampled_preamble_in_rx_filt = start_idx_oversampled_preamble_in_rx_filt + (num_preamble_symbols * OSR) - 1;
        if end_idx_oversampled_preamble_in_rx_filt <= length(rx_signal_filtered_complex) && start_idx_oversampled_preamble_in_rx_filt > 0
            rx_preamble_part_oversampled = rx_signal_filtered_complex(start_idx_oversampled_preamble_in_rx_filt : end_idx_oversampled_preamble_in_rx_filt);
            downsampling_start_index_for_cfo_est = 1 + best_known_sample_offset;
            rx_preamble_syms_for_cfo = [];
            if downsampling_start_index_for_cfo_est <= length(rx_preamble_part_oversampled)
                rx_preamble_syms_for_cfo = rx_preamble_part_oversampled(downsampling_start_index_for_cfo_est : OSR : end);
                if length(rx_preamble_syms_for_cfo) > num_preamble_symbols
                    rx_preamble_syms_for_cfo = rx_preamble_syms_for_cfo(1:num_preamble_symbols);
                end
            end
            if length(rx_preamble_syms_for_cfo) == num_preamble_symbols
                sum_corr_cfo = 0;
                for k_cfo = 0:(num_preamble_symbols/2 - 1)
                    idx1 = 2*k_cfo + 1; idx2 = 2*k_cfo + 2;
                    sum_corr_cfo = sum_corr_cfo + rx_preamble_syms_for_cfo(idx2) * conj(rx_preamble_syms_for_cfo(idx1));
                end
                delta_f_estimated = angle(sum_corr_cfo) / (2 * pi * Ts);
            end
        end
    end
    estimated_cfo_values(j_idx) = delta_f_estimated; % << STORE ESTIMATED CFO
    if mod(j_idx, round(length(SNR_dB_range)/5)) == 1 || j_idx == length(SNR_dB_range) % Print periodically
        fprintf('  Actual CFO: %.2f Hz, Estimated CFO: %.2f Hz\n', delta_f_actual_setting, delta_f_estimated);
    end

    % --- CFO CORRECTION ---
    if ~enable_cfo_correction % If correction disabled, force estimate to 0 for correction math
        delta_f_to_correct = 0;
    else
        delta_f_to_correct = delta_f_estimated;
    end
    num_rx_filt_samples = length(rx_signal_filtered_complex);
    t_osr_rx_filt = (0:num_rx_filt_samples-1) * (Ts/OSR); 
    rx_signal_cfo_corrected = rx_signal_filtered_complex .* exp(-1j * 2 * pi * delta_f_to_correct * t_osr_rx_filt);
    
    % ==== Downsample, Isolate Payload, Demodulate, BER & EVM ====
    % (Logic from your previous working script - ensure this part is robust)
    if length(rx_signal_cfo_corrected) > total_filter_delay_samples
        rx_symbols_delayed_complex_AC = rx_signal_cfo_corrected(1 + total_filter_delay_samples : end);
    else, rx_symbols_delayed_complex_AC = complex([]); end
    downsampling_start_index_data = 1 + best_known_sample_offset; 
    if downsampling_start_index_data <= length(rx_symbols_delayed_complex_AC)
        rx_symbols_downsampled_complex_AC = rx_symbols_delayed_complex_AC(downsampling_start_index_data : OSR : end);
    else, rx_symbols_downsampled_complex_AC = complex([]); end

    rx_payload_symbols_final = complex([]); 
    tx_payload_symbols_for_evm_AC = complex([]);
    original_payload_bits_for_ber_AC = [];
    start_payload_idx_in_downsampled = num_preamble_symbols + num_guard_symbols + 1;

    if length(rx_symbols_downsampled_complex_AC) >= (start_payload_idx_in_downsampled -1 + num_payload_symbols)
         potential_rx_payload_AC = rx_symbols_downsampled_complex_AC(start_payload_idx_in_downsampled : start_payload_idx_in_downsampled + num_payload_symbols - 1);
    elseif length(rx_symbols_downsampled_complex_AC) >= start_payload_idx_in_downsampled 
        potential_rx_payload_AC = rx_symbols_downsampled_complex_AC(start_payload_idx_in_downsampled : end);
    else, potential_rx_payload_AC = complex([]); end
    num_rx_payload_to_compare_AC = min(length(potential_rx_payload_AC), num_payload_symbols);

    if num_rx_payload_to_compare_AC > 0
        rx_payload_symbols_final = potential_rx_payload_AC(1:num_rx_payload_to_compare_AC);
        tx_payload_symbols_for_evm_AC = payload_symbols_qpsk_complex(1:num_rx_payload_to_compare_AC);
        original_payload_bits_for_ber_AC = payload_data_bits(1 : num_rx_payload_to_compare_AC * BitsPerSymbol_QPSK);
        
        if ~exist('QPSK_dem','file'), error('QPSK_dem.m not found.'); end
        rx_payload_bits = QPSK_dem(real(rx_payload_symbols_final), imag(rx_payload_symbols_final));
        len_rx_bits = length(rx_payload_bits);
        if len_rx_bits == length(original_payload_bits_for_ber_AC) && len_rx_bits > 0
            num_errors = sum(rx_payload_bits ~= original_payload_bits_for_ber_AC);
            BER(j_idx) = num_errors / len_rx_bits;
        else, BER(j_idx) = 1; end

        P_ref = mean(abs(tx_payload_symbols_for_evm_AC).^2); 
        if P_ref < 1e-9, P_ref = 1; end
        error_vector = rx_payload_symbols_final - tx_payload_symbols_for_evm_AC;
        evm = sqrt(mean(abs(error_vector).^2) / P_ref);
        evm_percent(j_idx) = evm * 100;
    else
        BER(j_idx) = 1; evm_percent(j_idx) = 100;
    end
    fprintf('  BER: %.2e, EVM: %.2f%%\n', BER(j_idx), evm_percent(j_idx));
end % End of SNR loop

%% Plot BER
figure;
semilogy(SNR_dB_range, BER, 'bo-', 'LineWidth', 1.5); grid on;
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate (BER)');
title_str = sprintf('QPSK BER (Actual CFO: %.0f Hz, Correction: %s, Offset: %d)', ...
                    delta_f_actual_setting, mat2str(enable_cfo_correction), best_known_sample_offset);
title(title_str); ylim_min_ber = max(1e-6, min(BER(BER>0))/2); if isempty(ylim_min_ber), ylim_min_ber = 1e-6; end
ylim([ylim_min_ber 1]); legend('Simulated BER','Location','southwest');

%% Plot EVM
figure;
plot(SNR_dB_range, evm_percent, 'ro-', 'LineWidth', 1.5); grid on;
xlabel('Eb/No (dB)'); ylabel('EVM (%)');
title_str_evm = sprintf('QPSK EVM (Actual CFO: %.0f Hz, Correction: %s, Offset: %d)', ...
                    delta_f_actual_setting, mat2str(enable_cfo_correction), best_known_sample_offset);
title(title_str_evm); legend('Simulated EVM','Location','northwest');

%% << NEW PLOT: Estimated CFO vs. SNR >>
figure;
plot(SNR_dB_range, estimated_cfo_values / 1000, 'mx-', 'LineWidth', 1.5, 'MarkerSize', 8); % Plot in kHz
hold on;
line(SNR_dB_range, ones(size(SNR_dB_range)) * delta_f_actual_setting / 1000, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('Eb/No (dB)');
ylabel('Estimated CFO (kHz)');
title_str_est_cfo = sprintf('Estimated CFO vs. SNR (Actual CFO: %.0f Hz = %.1f kHz)', ...
                            delta_f_actual_setting, delta_f_actual_setting/1000);
title(title_str_est_cfo);
legend('Estimated CFO', 'Actual CFO', 'Location', 'best');

%% Plot constellation of last SNR point 
if exist('rx_payload_symbols_final', 'var') && ~isempty(rx_payload_symbols_final)
    constDiagLast = comm.ConstellationDiagram(...
        'Title', ['Received Payload Constellation (Last Eb/No = ', num2str(current_EbNo_dB), ' dB)'], ...
        'XLimits', [-1.5, 1.5], 'YLimits', [-1.5, 1.5], ... 
        'ReferenceConstellation', refConstellation_qpsk_data, 'ShowLegend', true);
    step(constDiagLast, rx_payload_symbols_final'); 
else
    disp('No final payload symbols to plot constellation for the last SNR point.');
end
fprintf('\nMain BER simulation finished.\n');
