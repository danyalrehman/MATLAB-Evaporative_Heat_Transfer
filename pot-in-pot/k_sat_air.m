function k_mix = k_sat_air(T_Kelvin)

% This calculates the thermal conductivity of saturated air at a given
% temperature 

% Input is temperature celsius of the saturated air (air + water gas
% mixture) IN KELVIN

% Correlation Parmaters: https://pdfs.semanticscholar.org/cec4/877708079f2c645c5063508403e2f17163eb.pdf


% Convert to Celsius for the correlation

T_celsius = T_Kelvin - 273.15; % Celsius 

SK_0 = 2.40073953E-2;
SK_1 = 7.278410162E-5;
% SK_2 = 1.788037411E-2% likely incorrectly reported in the paper  
SK_3 = -1.351703529E-9;
SK_4 = -3.322412767E-11; 


k_mix = SK_0 + SK_1*T_celsius + SK_3.*T_celsius^3 + SK_4.*T_celsius^4; % W/mK


end