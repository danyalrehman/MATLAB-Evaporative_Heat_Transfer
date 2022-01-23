function Sky_temp = T_sky(T_inf, e_clear, C_a)


% T_inf in Celsius

T_inf = T_inf - 273.15; 

Sky_temp = (C_a^(0.25))*(e_clear^(0.25))*T_inf; % Celsius

% Convert to Kelvin

Sky_temp = Sky_temp + 273.15; % Kelvin



end