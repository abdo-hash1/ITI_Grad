function[noisy_signal]= adding_awgn_noise(signal,SNR)

 signal_power=mean(abs(signal(:)).^2);
% fprintf('\n\nSignal Power before adding awgn =%.10f ...\n', signal_power);
  SNR_linear = 10.^(SNR/10); 
  
noise_power= signal_power ./ SNR_linear;

if(isreal(signal))
    awgn_noise = sqrt(noise_power) * randn(size(signal));
else
    awgn_noise = sqrt(noise_power/2) * ( randn(size(signal)) + 1j * randn(size(signal)) );
        % awgn_noise = sqrt(noise_power/2) * (noise + 1j * noise) ;

end

noisy_signal = signal + awgn_noise;
% fprintf('\n\nSignal Power after adding AWGN (noisy_signal = signal + awgn_noise) =%.10f ...\n', mean(abs(noisy_signal(:)).^2));

end
