function T_dew_point = T_dp(T_inf, RH)

% T_inf is the ambient temperature in Kelvin
% RH is the relative humidity as a fraction 

% dew point calculation from w
% pw=(P*w)/(0.621945+w); % water vapor partial pressure in kPa

% Partial Pressure of Water
Pw = RH*Saturation_pressure(T_inf); % Pa 

% Convert to kPa for this correlation 

Pw = Pw/1000; % kPa

alpha=log(Pw);

T_dew_point =6.54 + 14.526*alpha + 0.7389*(alpha^2) + 0.09486*(alpha^3) + ...
    0.4569*(Pw^0.1984); % valid for Tdp between 0 C and 93 C

% Convert to Kelvin 



end