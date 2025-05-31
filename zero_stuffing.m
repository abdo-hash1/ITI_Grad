function [upsampled] = zero_stuffing(data,ratio)

N = length(data) ; 

upsampled = zeros(1,N*ratio);

upsampled(1:ratio:end) = data; 



end