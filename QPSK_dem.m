function demod_data = QPSK_dem(I,Q)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

demod_data = [] ; 
I(I<0) = 0;
Q(Q<0) = 0;
I(I>0) = 1;
Q(Q>0) = 1;
for i = 1:length(I)
    demod_data = [demod_data I(i) Q(i)] ;
end
demod_data(demod_data >0) = 1 ;
end