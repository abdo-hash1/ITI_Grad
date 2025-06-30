function [z_acc] = CORDIC_Vector(vector)

    CORDIC_ANGLES_deg = [
        45.0, 26.56505118, 14.03624347, 7.12501635, 3.57633437, 1.78991063, ...
        0.89517371, 0.44761417, 0.22381053, 0.11190568, 0.05595293, 0.02797650, ...
        0.01398826, 0.00699413, 0.00349706, 0.00174853, 0.00087427, 0.00043713, ...
        0.00021857, 0.00010928, 0.00005464, 0.00002732, 0.00001366, 0.00000683
    ];

x = real(vector);
y = imag(vector);
z_acc = 0 ;

for i = 1:length(CORDIC_ANGLES_deg)
    x_temp = x; y_temp = y;

    if y >= 0 
        d = -1;
    else 
        d = 1 ; 
    end 
    
    x = x_temp - d * y_temp* (2^(-(i-1))) ;
    y = y_temp + d * x_temp *  (2^(-(i-1))) ;
    
    z_acc = z_acc - d * CORDIC_ANGLES_deg(i);

    
end


end

