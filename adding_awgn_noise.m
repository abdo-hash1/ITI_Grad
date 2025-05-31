function[noisy_signal]= adding_awgn_noise(signal,SNR)

signal_power=mean(abs(signal).^2);
noise_power=signal_power/(10^(SNR/10));

if(isreal(signal))
    awgn_noise = sqrt(noise_power) * randn(size(signal));
else
    awgn_noise = sqrt(noise_power/2)*(randn(size(signal)) + 1j * randn(size(signal)));
end

noisy_signal = signal + awgn_noise;

end
