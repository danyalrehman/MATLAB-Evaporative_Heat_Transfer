function q_irr_flux = q_abs(t)

% t is inputted as the fraction of the day, as calculated in Balances 
% Returns the solar irradiance as a function of time 
data = xlsread('solar_irradiation_formatted.xls');

t_data = data(2:end,3)./(3600*24); % seconds 

q_flux_data = data(2:end,2); % W/m2

q_irr_flux = interp1(t_data, q_flux_data, t); % W/m2 

end