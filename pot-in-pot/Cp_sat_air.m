function Cp_mix = Cp_sat_air(T_Kelvin)

% Input is temperature celsius of the saturated air (air + water gas
% mixture) IN KELVIN

% This function calculates the specific heat capacity of the gaseous
% mixture of air and water at saturated conditions (100% RH) 

% Correlation Parmaters: https://pdfs.semanticscholar.org/cec4/877708079f2c645c5063508403e2f17163eb.pdf
SC_0 = 1.00457;
SC_1 = 2.05063E-3;
SC_2 = -1.63153E-4;
SC_3 = 6.2123003E-6;
SC_4 = -8.83047888E-8;
SC_5 = 5.071307038E-10;

% Convert to Celsius for the correlation

T_celsius = T_Kelvin - 273.15; % Celsius 

Cp_mix = SC_0 + SC_1.*T_celsius + SC_2.*T_celsius.^2 + SC_3.*T_celsius.^3 + ...
    SC_4.*T_celsius.^4 + SC_5.*T_celsius.^5; % kJ/kgK

Cp_mix = Cp_mix*1000; % J/kgK


end