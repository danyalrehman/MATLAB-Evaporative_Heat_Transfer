function k_water = k_w(T_Kelvin)

% T_kelvin is temperature in Kelvin

% This function calculates the thermal conductivity of pure water as a function
% of temperature at atmospheric pressure (1 bar) 

% Source : https://www.nist.gov/sites/default/files/documents/srd/jpcrd38200921p.pdf

% Correlation requires temperature to be divided by 300K 
T = T_Kelvin./300;

k_water = 0.80201.*T.^(-0.32) + -0.25992.*T.^(-5.7) + ...
    0.10024.*T.^(-12.0) + -0.0032005.*T.^(-15.0); % W/mK 



end