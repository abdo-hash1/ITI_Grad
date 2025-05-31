clear;clc
% *******************************************
% ITI DSP Graduation Project
% *******************************************

% ==== System PARAMETERS ====
Rs = 1.92e6 ;
Ts = 1/Rs;     
alpha = 0.5;       
N_span = 6;         % Filter span (symbols)
OSR = 4;             % Upsampling factor
N= 144*1000;
SNR = 0:20;
BER = zeros(size(SNR));


%% Design an RRC filter
[h_filter, t_filter] = designPulseShapingFilter(Ts, alpha, N_span, OSR, 'RRC');

%%
% Plot the RC filter
figure;
stem(t_filter, h_filter, 'filled');
hold on;
plot(t_filter, h_filter, '-');
hold off;
title(['RRC Filter Impulse Response (alpha = ', num2str(alpha), ', OSR = ', num2str(OSR), ')']);
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
legend('Samples', 'Envelope', 'Location', 'best');
%%

for j = 1:length(SNR)

% ==== Generate random symbols (e.g., QPSK) ====

data = randi([0 1],1,N);
[I,Q] = QPSK(data);
% ==== Upsample symbols ====
I_upsampled = zero_stuffing(I,OSR);
Q_upsampled = zero_stuffing(Q,OSR);
% ==== Apply RRC Filter ====
I_filtered = conv(I_upsampled,h_filter,'same') ;
Q_filtered = conv(Q_upsampled,h_filter,'same') ;

rx_signal =  I_filtered + 1j* Q_filtered ;

noisy_signal = adding_awgn_noise(rx_signal,SNR(j));

% ==== Receiver ====
I_rx = conv(real(noisy_signal), h_filter, 'same');
Q_rx = conv(imag(noisy_signal), h_filter, 'same');

% ==== Downsample ===
I_down = I_rx(1:OSR:end);
Q_down = Q_rx(1:OSR:end);

% ==== Demodulate ====
rx_data = QPSK_dem(I_down, Q_down);

% ==== BER ====
num_errors = sum(rx_data ~= data);
BER(j) = num_errors / length(data);

% ==== EVM Calculation ====
tx_symbols = I + 1j*Q;
rx_symbols = I_down + 1j*Q_down;
rx_symbols = rx_symbols / sqrt(mean(abs(rx_symbols).^2)); 
tx_symbols = tx_symbols / sqrt(mean(abs(tx_symbols).^2));

evm = sqrt(mean(abs(rx_symbols - tx_symbols).^2))/ sqrt(mean(abs(tx_symbols).^2));
evm_percent(j) = evm * 100;
end

%%
semilogy(SNR, BER(i,:), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
title(['SNR vs. BER for Span ', num2str(N_span(i)) ]);
