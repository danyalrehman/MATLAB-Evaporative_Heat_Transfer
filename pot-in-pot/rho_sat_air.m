function rho_mix = rho_sat_air(T_Kelvin)

% This calculates the density of saturated air at a given
% temperature 

% Input is temperature celsius of the saturated air (air + water gas
% mixture) IN KELVIN

% Correlation Parmaters: https://pdfs.semanticscholar.org/cec4/877708079f2c645c5063508403e2f17163eb.pdf


% Convert to Celsius for the correlation

T_celsius = T_Kelvin - 273.15; % Celsius 


SD_0 = 1.293393662; 
SD_1 = -5.538444326E-3; 
SD_2 = 3.860201577E-5;
SD_3 = -5.2536065E-7; 

rho_mix = SD_0 + SD_1.*T_celsius + SD_2.*T_celsius.^2 + ...
    SD_3.*T_celsius.^3; % kg/m3 


end