function [h,t] = RRC(Ts,alpha,N_span,OSR)
   
    N_taps = N_span * OSR + 1 ;
    step = Ts / OSR;            % Time between samples
    
    % Time vector (centered at 0)
    t = (-(N_taps-1)/2 : (N_taps-1)/2) * step;
    
    % Initialize impulse response
    h = zeros(size(t));
    
    % Handle singularities
    for i = 1:length(t)
        ti = t(i);
        
        % Case 1: Center tap (t = 0)
        if ti == 0
            h(i) = (1 + alpha * (4/pi - 1)) / sqrt(Ts);
            
        % Case 2: Singularity at t = ±Ts/(4α)
        elseif abs(abs(ti) - Ts/(4*alpha)) < 1e-10  % Floating-point tolerance
            h(i) = (alpha / sqrt(2*Ts)) * ...
                   ((1 + 2/pi) * sin(pi/(4*alpha))) + ...
                   ((1 - 2/pi) * cos(pi/(4*alpha)));
            
        % Case 3: General case
        else
            term1 = sin(pi * ti/Ts * (1 - alpha));
            term2 = 4 * alpha * ti/Ts * cos(pi * ti/Ts * (1 + alpha));
            numerator = term1 + term2;
            
            denominator = pi * ti/Ts * (1 - (4 * alpha * ti/Ts)^2);
            
            h(i) = numerator / (denominator * sqrt(Ts));
        end
    end
    
    % Normalize to peak amplitude = 1
    h = h / max(abs(h));
end
