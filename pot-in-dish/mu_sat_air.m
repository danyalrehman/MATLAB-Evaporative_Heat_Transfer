function mu_mix = mu_sat_air(T_Kelvin)

% This calculates the viscosity of saturated air at a given
% temperature 

% Input is temperature celsius of the saturated air (air + water gas
% mixture) IN KELVIN

% Correlation Parmaters: https://pdfs.semanticscholar.org/cec4/877708079f2c645c5063508403e2f17163eb.pdf


% Convert to Celsius for the correlation

T_celsius = T_Kelvin - 273.15; % Celsius 

SV_0 = 1.715747771E-5;
SV_1 = 4.722402075E-8;
SV_2 = -3.663027156E-10;
SV_3 = 1.873236686E-12;
SV_4 = -8.050218737E-14;

mu_mix = SV_0 + SV_1*T_celsius + SV_2*T_celsius^2 + SV_3*T_celsius^3 + ...
    SV_4*T_celsius^4; % [Ns/m2] = [kg*m/s^2 * s/m2] --> [kg/ms] = [Pa*s] 

end