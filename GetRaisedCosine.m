function [h, t] = GetRaisedCosine(Tsymbol,alpha, N, Oversampling, Root) 
%Tsymbol = 1; %% Default = 1 second 

Step = 1/Oversampling; 

if(mod(N,2) == 1) 
    t = -floor(N/2)*Step : Step : floor(N/2)*Step; %% for odd N 
else 
    t = (-N/2 +0.5)*Step : Step : (N/2 -0.5)*Step;  
end 
t = t + 1e-8; % we add this small offset to avoid divide by zero errors 
%Evaluate Raised Cosine Pulse Shape at t 
h_RC = sinc(t/Tsymbol) .* ( cos(alpha*pi*t/Tsymbol)./(1 - (2*alpha*t/Tsymbol).^2));

% Evaluate Root Raised Cosine Pulse Shape at t 
RRC = (4*alpha ./ (pi*sqrt(Tsymbol)*(1-(4*alpha*t/Tsymbol).^2)) ) .* ... 
( cos((1+alpha)*pi*t/Tsymbol) + (Tsymbol./(4*alpha*t)).*sin((1- ...
alpha)*pi*t/Tsymbol) ); 

h_RRC = RRC/max(RRC); 
if(Root == 1)
    h = h_RRC; 
else 
    h = h_RC; 
end
end