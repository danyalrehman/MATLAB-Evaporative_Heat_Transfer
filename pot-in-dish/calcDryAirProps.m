function [dry_air_props] = calcDryAirProps(T_ambient_data)

% Assumed 0% relative humidity (Dry air)

% T_ambient_data is in Celsius and assumes a column vector input

T_ambient_data = transpose(T_ambient_data);

% Source: https://pdfs.semanticscholar.org/cec4/877708079f2c645c5063508403e2f17163eb.pdf

T = T_ambient_data + 273.15; 

% All correlations are correct with temperature input in Kelvin, 
%   despite the paper erroniously stating celsius 

% Viscosity Correlation 
MA_0 = -9.8601E-1;
MA_1 = 9.080125E-2;
MA_2 = -1.17635575E-4;
MA_3 = 1.2349703E-7;
MA_4 = -5.7971299E-11;

mu_a = MA_0 + MA_1.*T + MA_2.*T.^2 + MA_3.*T.^3 + MA_4.*T.^4; % [Pa*s * 10^(-6)]
mu_a = mu_a*10^(-6); % [Pa*s]


% Thermal Conductivity
KA_0 = -2.276501E-3;
KA_1 = 1.2598485E-4;
KA_2 = -1.4815235E-7;
KA_3 = 1.73550646E-10;
KA_4 = -1.066657E-13;
KA_5 = 2.47663035E-17;

k_a = KA_0 + KA_1.*T + KA_2.*T.^2 + KA_3.*T.^3 + KA_4.*T.^4 + ...
    KA_5.*T.^5; % [Pa*s * 10^(-6)]

% Heat Capacity 
CA_0 = 1.03409;
CA_1 = -0.284887E-3;
CA_2 = 0.7816818E-6;
CA_3 = -0.4970786E-9;
CA_4 = 0.1077024E-12;

Cp_a = CA_0 + CA_1.*T + CA_2.*T.^2 + CA_3.*T.^3 + CA_4.*T.^4; % 
Cp_a = Cp_a*1000; % J/kgK

% Density 

P_O = 101325; % Pa (atmospheric pressure)
R = 8.3134; % J/molK
M_a = 28.963; % kg/kmol

rho_a = 1E-3.*(M_a.*P_O)./(R.*T); % kg/m3 


dry_air_props = [mu_a; k_a; Cp_a; rho_a];


end