function q_irr_flux = q_abs(t)

% t is inputted as the fraction of the day, as calculated in Balances 
% Returns the solar irradiance as a function of time 
%data = readmatrix('solar_irradiation_formatted.xls');
loadedfile = load("solar_data.mat", "-mat");

data = loadedfile.data;
% Preload the data to save compute time 

t_data = data(1:end,3)./(3600*24); % seconds 

q_flux_data = data(1:end,2); % W/m2

q_irr_flux = interp1(t_data, q_flux_data, t); % W/m2 
end