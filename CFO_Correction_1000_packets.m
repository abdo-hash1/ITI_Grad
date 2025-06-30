% --- QPSK Transceiver with Coarse CFO ---

clear; clc; 

% --- 1. System Parameters ---
Rs = 1.92e6 ;           % Symbol rate
Ts = 1/Rs;
alpha = 0.5;            % Roll-off factor for RRC filter
N_span = 6;             % Filter span (symbols)
OSR = 4;                % Oversampling factor
N_payload_bits = 2240;  % Fixed payload size in bits
total_payload_data_tx_bits = randi([0 1], 1000, N_payload_bits);

SNR_dB_range = 0:1:60;  

% ---- CFO PARAMETER (using PPM) ----

% Choose a nominal RF carrier frequency for your system (in Hz)
fc_hz = 2.4e9; 

% Define your desired ppm value
ppm_offset = 25;

% Calculate the absolute frequency offset in Hz
max_offset_hz = (ppm_offset * 1e-6) * fc_hz;
step_hz = 1;    

% Create the vector of frequency offsets
cfo_range_hz = -max_offset_hz : step_hz : max_offset_hz;

% =======================================================================
% --- NEW: COSTAS PLL PARAMETERS ---
% =======================================================================
pll_loop_bandwidth = 0.02;  % Normalized loop bandwidth (adjust as needed)
pll_damping_factor = 0.707; % Damping factor (0.707 is a good starting point)

% Calculate Proportional (Kp) and Integral (Ki) gains for the PI loop filter
% This is a standard formula for converting bandwidth/damping to gains
theta = pll_loop_bandwidth / ( (pll_damping_factor + 1/(4*pll_damping_factor)) * 1 );
Kp = (4 * pll_damping_factor * theta) / (1 + 2*pll_damping_factor*theta + theta^2);
Ki = (4 * theta^2) / (1 + 2*pll_damping_factor*theta + theta^2);



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
    for l = 1:1000
    random_cfo_value = cfo_range_hz(randi(length(cfo_range_hz)));
    current_SNR = SNR_dB_range(j_idx); 
    fprintf('\n\nSimulating for SNR = %.1f dB and CFO =%.1f ...\n', current_SNR,random_cfo_value);
    % ==== Transmit Path ==== 
    % Transmitted Bits
    payload_data_tx_bits = total_payload_data_tx_bits(l,:);
    BitsPerSymbol = 2;
    % I/Q Modulation
    [I_payload, Q_payload] = QPSK(payload_data_tx_bits);
    % Adding Preamble to Payload 
    payload_symbols_qpsk_complex = I_payload + 1j*Q_payload;
    num_payload_symbols = length(payload_symbols_qpsk_complex);
    total_symbols_complex_tx = [preamble_symbols_complex,payload_symbols_qpsk_complex];
    % Upsampling

    upsampled_total_symbols_complex = upsample(total_symbols_complex_tx, OSR);
    % Pulse Shaping
    tx_signal_shaped = conv(upsampled_total_symbols_complex, h_filter);

    % ==== Channel Model ==== 
    % Adding CFO
    num_tx_samples = length(tx_signal_shaped);
    t_osr_tx = (0:num_tx_samples-1) * (Ts/OSR);
    tx_signal_with_cfo = tx_signal_shaped .* exp(1j * 2 * pi * random_cfo_value * t_osr_tx);
    % Adding AWGN Noise
    snr_for_awgn_func_dB = current_SNR + 10*log10(BitsPerSymbol / OSR);
    % noisy_signal = awgn(tx_signal_with_cfo, snr_for_awgn_func_dB, 'measured');
    noisy_signal = adding_awgn_noise(tx_signal_with_cfo,snr_for_awgn_func_dB) ;
    % ==== RECEIVER ====
   rx_signal_filtered_complex = conv(noisy_signal, h_filter);

    % Downsampling & Truncation
    total_delay = length(h_filter) - 1;
    rx_downsampled_long = rx_signal_filtered_complex(total_delay + 1 : OSR : end);
    num_tx_symbols = length(total_symbols_complex_tx);
    rx_downsampled = rx_downsampled_long(1:num_tx_symbols);
    
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
    fprintf('Packet %.2f |  Actual CFO: %.2f Hz | Manual Estimate: %.2f Hz | Difference: %.2f Hz  \n', l ,random_cfo_value, delta_f_hat, abs(random_cfo_value-delta_f_hat));
    
    % Coarse CFO Correction (applied to the whole packet)
    t_full_packet = (0:length(rx_downsampled)-1).' * Ts;
    cfo_corrector_full = exp(-1j * 2 * pi * delta_f_hat * t_full_packet);
    rx_packet_coarse_corrected = rx_downsampled(:) .* cfo_corrector_full(:);
    
    
    %  Extract the Payload for Demodulation
    payload_symbols_for_demod = rx_packet_coarse_corrected(N+1:end);
    num_payload_symbols = length(payload_symbols_for_demod);


    % =======================================================================
    % --- NEW: COSTAS PLL for Fine Frequency and Phase Tracking ---
    % =======================================================================
    % Initialize PLL state variables
    pll_phase = 0;
    integrator_out = 0;
    payload_symbols_after_pll = zeros(num_payload_symbols, 1);

    for n = 1:num_payload_symbols
        % 1. Get input symbol
        in_sym = payload_symbols_for_demod(n);
        
        % 2. Apply phase correction from NCO
        corrected_sym = in_sym * exp(-1j * pll_phase);
        
        % 3. Store the fully corrected symbol
        payload_symbols_after_pll(n) = corrected_sym;
        
        % 4. Phase Error Detector (for QPSK)
        %    This calculates how far the corrected symbol is from the ideal axes.
        phase_error = sign(real(corrected_sym)) * imag(corrected_sym) - ...
                      sign(imag(corrected_sym)) * real(corrected_sym);
        
        % 5. Loop Filter (Proportional-plus-Integrator)
        integrator_out = integrator_out + Ki * phase_error;
        proportional_out = Kp * phase_error;
        
        % 6. Update NCO phase for the *next* symbol
        pll_phase = pll_phase + proportional_out + integrator_out;
    end
    % =======================================================================
    
    
    % --- Demodulation ---
    rx_payload_bits = QPSK_dem(real(payload_symbols_after_pll), imag(payload_symbols_after_pll));
    % ==== BER Calculation on PAYLOAD ====
    num_errors = sum(rx_payload_bits ~= payload_data_tx_bits);
    BER_simulated_packet(l) = num_errors / N_payload_bits;
    end

    BER_simulated(j_idx) = sum(BER_simulated_packet(l)) / 1000;

    % fprintf('  BER: %.2e\n', BER_simulated(j_idx));

    % ==== Theoretical BER ====
    ebno_linear = 10^(current_SNR/10);
    BER_theoretical(j_idx) = qfunc(sqrt(2 * ebno_linear));
end % End of SNR loop

%% Plotting
% BER Plot
figure;
semilogy(SNR_dB_range, BER_simulated, 'bo-', 'LineWidth', 1.5);
hold on; semilogy(SNR_dB_range, BER_theoretical, 'r--', 'LineWidth', 1.5);
hold off; grid on; xlabel('Eb/No (dB)'); ylabel('BER');
title(sprintf('QPSK BER with CFO Correction'));
legend('Simulated BER', 'Theoretical BER', 'Location','southwest');
ylim_min_ber = max(1e-6, min(BER_simulated(BER_simulated>0))/2); if isempty(ylim_min_ber), ylim_min_ber = 1e-6; end; ylim([ylim_min_ber 1]);
%%
% Constellation Plot of last SNR point
if exist('rx_payload_symbols_before_pll', 'var') && ~isempty(rx_payload_symbols_before_pll) && ...
   exist('rx_payload_symbols_after_pll', 'var') && ~isempty(rx_payload_symbols_after_pll)
    figure;
    subplot(1,2,1);
    plot(rx_payload_symbols_before_pll(1:min(end,2000)),'.');
    title('Payload Symbols BEFORE Fine PLL'); axis equal; grid on; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    
    subplot(1,2,2);
    plot(rx_payload_symbols_after_pll(1:min(end,2000)),'.');
    title('Payload Symbols AFTER Fine PLL'); axis equal; grid on; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    sgtitle(sprintf('Constellations at Last Eb/No = %.1f dB', current_SNR));
else
    disp('No final payload symbols to plot constellation.');
end
fprintf('\nSimulation finished.\n');

% Estimate over different cfo ( min to max)
