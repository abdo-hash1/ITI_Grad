function [h, t] = designPulseShapingFilter(Ts, alpha, N_span, OSR, filterType)
    % --- Filter Parameters ---
    Tsample = Ts / OSR;         
    Ntaps = N_span * OSR + 1;   
   
    % Create time vector centered around 0
    % It spans from -N_span*Ts/2 to +N_span*Ts/2
    t = (-(Ntaps-1)/2 : (Ntaps-1)/2) * Tsample;

    h = zeros(1, Ntaps);       
    tol = 1e-12;                

    xt = t / Ts; 

    % --- Generate Filter Coefficients ---
    if strcmpi(filterType, 'RC')
        for i = 1:Ntaps
            x_curr = xt(i); 
            if abs(alpha) < tol 
                h(i) = sinc(x_curr);
            else
                if abs(x_curr) < tol 
                    h(i) = 1.0;
                elseif abs(abs(x_curr) - 1/(2*alpha)) < tol 
                    h(i) = (alpha/2) * sin(pi/(2*alpha));
                else
                    % Standard RC formula
                    numerator_rc = cos(pi * alpha * x_curr);
                    denominator_rc = (1 - (2 * alpha * x_curr)^2);
                    if abs(denominator_rc) < tol 
                        h(i) = 0; 
                    else
                        h(i) = sinc(x_curr) * (numerator_rc / denominator_rc);
                    end
                end
            end
        end

    elseif strcmpi(filterType, 'RRC')
        for i = 1:Ntaps
            x_curr = xt(i); 

            if abs(alpha) < tol 
                h(i) = sinc(x_curr); 
            else
                if abs(x_curr) < tol
                    h(i) = (1 - alpha + 4*alpha/pi);
                elseif abs(abs(x_curr) - 1/(4*alpha)) < tol 
                    h(i) = (alpha/sqrt(2)) * ( (1 + 2/pi)*sin(pi/(4*alpha)) + ...
                                               (1 - 2/pi)*cos(pi/(4*alpha)) );
                else
                    % Standard RRC formula
                    numerator_rrc = sin(pi*x_curr*(1-alpha)) + 4*alpha*x_curr .* cos(pi*x_curr*(1+alpha));
                    denominator_rrc = pi*x_curr .* (1 - (4*alpha*x_curr).^2);
                    
                    if abs(denominator_rrc) < tol 
                        h(i) = 0; 
                    else
                        h(i) = numerator_rrc / denominator_rrc;
                    end
                end
            end
        end
        
        if any(h) 
           h = h / sqrt(sum(h.^2));
        end
    end
end

