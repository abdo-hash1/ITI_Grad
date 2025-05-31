clear;clc
% *******************************************
% ITI DSP Graduation Project
% *******************************************

% ==== System PARAMETERS ====
Rs = 1.92e6 ;
Ts = 1/Rs;     
alpha = 0.5;       
N_sym = 6;         % Filter span (symbols)
L = 4;             % Upsampling factor
N= 144*1000;
SNR = 0:20;
BER = zeros(size(SNR));

% ==== RRC Filter ====
[rrc, t_filter] = RRC(Ts, alpha, N_sym, L);
%%
%Plot impulse response
stem(t_filter, rrc, 'filled', 'MarkerSize', 4);
title('RRC Impulse Response');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;
%%
%Plot data after modulating 

t = (0:length(I)-1)*Ts;

figure;
subplot(2,1,1);
plot(t, I, 'b');
title('In-Phase Component (I)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 10*Ts]);  

subplot(2,1,2);
plot(t, Q, 'r');
title('Quadrature Component (Q)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 10*Ts]);  

figure;
plot(I, Q, '.');
axis equal;
grid on;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
title('QPSK Constellation Diagram');
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
%% 
% ===== Loop over the whole system ====
for idx = 1:length(SNR)

% ==== Generate random symbols (e.g., QPSK) ====

data = randi([0 1],1,N);
[I,Q] = QPSK(data);

% ==== Upsample symbols ====
I_upsampled = zero_stuffing(I,4);
Q_upsampled = zero_stuffing(Q,4);

% ==== Apply RRC Filter ====
I_filtered = conv(I_upsampled,rrc,'same') ;
Q_filtered = conv(Q_upsampled,rrc,'same') ;

rx_signal =  I_filtered + 1j* Q_filtered ;

noisy_signal = adding_awgn_noise(rx_signal,SNR(idx));

% ==== Receiver ====
I_rx = conv(real(noisy_signal), rrc, 'same');
Q_rx = conv(imag(noisy_signal), rrc, 'same');

% ==== Downsample ===
I_down = I_rx(1:L:end);
Q_down = Q_rx(1:L:end);

% ==== Demodulate ====
rx_data = QPSK_dem(I_down, Q_down);

% ==== BER ====
num_errors = sum(rx_data ~= data);
BER(idx) = num_errors / length(data);
end
%%
% EVM Calculation
tx_symbols = I + 1j*Q;
rx_symbols = I_down + 1j*Q_down;
rx_symbols = rx_symbols / sqrt(mean(abs(rx_symbols).^2)); 
tx_symbols = tx_symbols / sqrt(mean(abs(tx_symbols).^2));

evm = sqrt(mean(abs(rx_symbols - tx_symbols).^2))/ sqrt(mean(abs(tx_symbols).^2));
evm_percent = evm * 100;
fprintf('EVM = %.2f%%\n', evm_percent);
%%
rx_symbols = I_down + 1j*Q_down;
rx_symbols = rx_symbols / sqrt(mean(abs(rx_symbols).^2));

figure;
plot(real(rx_symbols), imag(rx_symbols), 'o');
grid on;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Constellation Diagram');
axis equal;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);