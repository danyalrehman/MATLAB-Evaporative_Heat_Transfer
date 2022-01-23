function [Ambient_air_props] = calcAirProps(T_ambient_data, ambient_humidity_data)

% This function calculates the thermal properties of ambient air as a
% function of ambient temperature and humidity 

% Source: https://pdfs.semanticscholar.org/cec4/877708079f2c645c5063508403e2f17163eb.pdf


% NOTE: T_Ambient_data is in CELSIUS

dry_air_props = calcDryAirProps(T_ambient_data); % 4 x n matrix

% Convert to Kelvin

T_ambient_data = T_ambient_data + 273.15; % Kelvin

P_O = 101325; % Pa (atmospheric pressure)
R = 8.3134; % J/molK
M_a = 28.963; % kg/kmol
M_v = 18.02; % kg/kmol

% Obtain the dry air properties (must insert a [1 x n] vector into
% calcDryAirProps function) 

mu_a = transpose(dry_air_props(1, :)); % Pa*s
k_a = transpose(dry_air_props(2, :)); % W/mK
Cp_a = transpose(dry_air_props(3, :)); % J/kgK
rho_a = transpose(dry_air_props(4, :)); % kg/m3

% Calculate the vapor mole fraction of water 
x_v = ambient_humidity_data.*Saturation_pressure(T_ambient_data)./P_O;

% Assume compressibility factor (z) --> 1 since T_ambient < 40 C for all cases 

%%  Calculate the density of the air/water mixture [kg/m3]
rho_ambient = M_a.*(P_O./(R.*(T_ambient_data))).*(1 - x_v.*(1 - (M_v./M_a)))./1000;  
% Column vector of rho values 

%% Calculate the Viscosity of the air/water mixture [Pa*s]
% Calculate the vapor mole fraction of air 
x_a = 1 - x_v; 

% Calculate Viscosity of water as a function of temperature 
% Source: https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
A = 2.414*1E-5; % Pa*s
B = 247.8; % Kelvin
C = 140; % Kelvin 

mu_v = A.*10.^(B./(T_ambient_data - C));

% Calculate interaction parameters
% mu_v = 8.90*1E-4; % Pa*s
phi_av = ((1 + ((mu_a./mu_v).^(1/2)).*(M_v./M_a)).^2)./((8.*(1 + ...
    M_a./M_v)).^(1/2));
phi_va = (mu_v./mu_a).*(M_a./M_v).*phi_av;

% Calculate the viscosity of the air/water mixture
mu_ambient = ((x_a.*mu_a)./(x_a + x_v.*phi_av)) + ...
    ((x_v.*mu_v)./(x_v + x_a.*phi_va)); % [Pa*s]

%%  Calculate the thermal conductivity of the air/water mixture 
k_v = k_w(T_ambient_data);

k_ambient = ((x_a.*k_a)./(x_a + x_v.*phi_av)) + ...
    ((x_v.*k_v)./(x_v + x_a.*phi_va)); % W/mK 


%% Heat Capacity of the air/water mixture (ambient) 
M_m = M_a.*x_a + M_v.*x_v;

Cp_v = 4184; % J/kgK % Can be taken to be constant in our range (Koretsky, 2nd Edition)  
Cp_m = Cp_a.*x_a.*(M_a./M_m) + Cp_v.*x_v.*(M_v./M_m); % J/kgK 


%% Output the properties

Ambient_air_props = [mu_ambient, k_ambient, Cp_m, rho_ambient];
% Returns a n x 4 array of ambient air propreties, where n is the amount of
% time steps you have humidity and temperature data for 

end