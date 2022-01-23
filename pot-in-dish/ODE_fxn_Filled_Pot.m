function [derivs] = ODE_fxn_Filled_Pot(t, vars, parameters, other)


derivs = zeros(length(vars), 1);

%% Define Variables 

% All in Kelvin 

% Temperature 

T_air = vars(1); % inner chamber, tomatoes and air mixture
T_veges = vars(2);
T_2 = vars(3); % clay pot
T_3 = vars(4); % sand 
T_4 = vars(5); % dish (probably plastic)
T_5 = vars(6); % cloth 
V_water = vars(7); % m3 of water in device 

T_roof_ext = vars(8); % Kelvin, external roof surface temp
T_roof_int = vars(9); % Kelvin, intenral roof surface temp 

h_cap = vars(10); % [m]

t_data = other(:, 1); % seconds 
T_outside = other(:, 2) + 273.15; % Kelvin
outside_humidity = other(:, 3); 
air_props = other(:, 4:7); 

mu_inf = interp1(t_data, air_props(:,1), t);
k_inf = interp1(t_data, air_props(:,2), t);
Cp_inf = interp1(t_data, air_props(:,3), t);
rho_inf = interp1(t_data, air_props(:,4), t);
nu_inf = mu_inf/rho_inf; % kinematic viscosity of ambient air 
alpha_inf = k_inf/(rho_inf*Cp_inf); % thermal diffusivity of ambient air

T_ambient = interp1(t_data, T_outside, t);

RH_inf = interp1(t_data, outside_humidity, t);

% Volumes
V_vec = parameters(:,3); % Volume 

V_air = V_vec(1);
V_veges = V_vec(2);
% V_1 = V_air + V_veges; % m3 
V_2 = V_vec(3);
V_3 = V_vec(4);
V_4 = V_vec(5);
V_5 = V_vec(6);

% Unpack all the parameters 
Cp_vec = parameters(:, 2); % Cp

Cp_air = Cp_sat_air(T_air);
Cp_veges = Cp_vec(2);
Cp_2 = Cp_vec(3);
Cp_3 = Cp_vec(4); % Sand 
Cp_4 = Cp_vec(5); % dish 
Cp_5 = Cp_vec(6);

Cp_water = 4184; % J/kgK
x_water = V_water/V_3; % volume fraction of water (max at 0.4, min at 0) 
sand_porosity = parameters(5, 9);

% Weighted average Cp of sand
Cp_3 = Cp_water*x_water + sand_porosity*Cp_3;


outer_speed = parameters(2, 9); 
inner_speed = parameters(1, 9); 
speeds = [inner_speed, outer_speed];


% Densities 
rho_vec = parameters(:,1); % Density
rho_air = rho_sat_air(T_air);
rho_veges = rho_vec(2);
rho_2 = rho_vec(3);
rho_3 = rho_vec(4);
rho_4 = rho_vec(5);
rho_5 = rho_vec(6);

% Thermal Conductivity 
k_vec = parameters(:, 6); % W/mK
k_air = k_sat_air(T_air);
k_veges = k_vec(2);
k2 = k_vec(3);
k3 = k_vec(4);
k4 = k_vec(5);
k5 = k_vec(6);

% Radii and Heights 

r_vec = parameters(:, 4);

r1 = r_vec(2);
r2 = r_vec(3);
r3 = r_vec(4);
r4 = r_vec(5);
t_s = r_vec(6); % sand depth below pot and above dish 

H_vec = parameters(:, 5);

H_air = H_vec(1);
H_veges = H_vec(2);
H1 = (H_air + H_veges); % full height 
H2 = H_vec(3);
H3 = H_vec(4);
H4 = H_vec(5);
L5 = H_vec(6);

% Heat Transfer Coefficients 
ambient_vars = [T_ambient, RH_inf, mu_inf, k_inf, Cp_inf, rho_inf];

h_vec = ConvectionModel(speeds, vars, t, parameters, ambient_vars);
h1 = h_vec(1);
h2 = h_vec(2);
h3 = h_vec(3);
h4 = h_vec(4);
h5 = h_vec(5);
 

% Calculate T_sky 

% Calcualte Ambeint Dew Point Temperature
% Cloud cover factor  (Clark and Allen) 
N = 0; % Clear skies 
CC = 1 + 0.0224*N - 0.0035*N^2 + 0.00028*N^3; % Correlation



Tdp_amb = T_dp(T_ambient, RH_inf);
Tsky = T_sky(T_ambient, E_clear(Tdp_amb), CC);

emissivity_vec = parameters(:, 8); 
E_sand = emissivity_vec(1);
E_dish = emissivity_vec(2);
E_clay = emissivity_vec(3);
E_cloth = emissivity_vec(4);
E_roof = emissivity_vec(5);
E_wall = emissivity_vec(6);

sigma = 5.6703*10^(-8); % W/m^2K^4
  
% Surface Areas
SA_vec = parameters(:,7);
A_1_bottom = SA_vec(2); % Area of contact between bottom of surface 1 and contacting part of
                        % surface 2 
                 
SA_2 = SA_vec(3);                      
SA_3 = SA_vec(4);              
SA_4 = SA_vec(5);
SA_5 = SA_vec(6);

sand_thickness = r3 - r2; 
dish_thickness = r4 - r3;


%% Radiation Networks and Roof/Dwelling calculations 


z_roof = parameters(1, 10); % [m]; Height from ground to the roof 
house_width = z_roof; % [m]
house_length = house_width; % [m] 

u_roof = parameters(2, 10); % [m/s]
A_roof = house_width*house_length; % m2

Re_roof = rho_inf*u_roof*(sqrt(A_roof))/(mu_inf);
Pr_inf = nu_inf/alpha_inf;

A_wall = z_roof*house_width; % [m2]
A_g = house_width*house_length - SA_5; % Area of the ground [m2]


% Assume only forced convection
Nu_ext_roof = 0.664*sqrt(Re_roof)*(Pr_inf^(1/3));
% Ra_roof = g*B_outside*(T_roof_ext - Tsky)*(L_roof^3)/(nu_inf*alpha_inf);
h_conv_roof_ext = Nu_ext_roof*(k_inf)/(sqrt(A_roof)); % [W/m2K]


h_rad_roof_ext = E_roof*sigma*(T_roof_ext + Tsky)*(T_roof_ext^2 + Tsky^2);  

% Conduction in the roof 
L_roof = parameters(3, 10); % [m]  thickness of the roof

k_roof = parameters(6, 10); % [W/mK] 
rho_roof = parameters(4, 10); % kg/m3
Cp_roof = parameters(5, 10); % [J/kgK]

Q_cond_roof = k_roof*A_roof*(T_roof_ext - T_roof_int)/L_roof;
Q_rad_ext = h_rad_roof_ext*A_roof*(Tsky - T_roof_ext);
Q_conv_ext = h_conv_roof_ext*A_roof*(T_ambient - T_roof_ext);

V_roof = L_roof*A_roof; % [m3] 

% Interior

% Surface 0 represents the inside walls and roof of the dwelling 

E_0 = (E_roof + E_wall)/2;
E_s = (E_clay + E_cloth + E_sand + E_dish)/2; 

A_0 = A_roof + 4*A_wall; 
A_s = SA_4 + SA_5 + SA_2 + SA_3; 

% Define: 

F_s0 = A_g/(A_0 + A_g);
F_sg = 1 - F_s0;
F_0g = 1 - F_s0*(A_s/A_0);

R_0 = (1 - E_0)/(A_0*E_0);
R_0s = 1/(A_s*F_s0);
R_s = (1 - E_s)/(A_s*E_s);
R_0g = 1/(A_0*F_0g);
R_sg = 1/(A_s*F_sg); 

R_eff_rad = R_0 + R_s + inv(inv(R_0s) + inv(R_0g + R_sg));

Q_rad_int = sigma*(T_roof_int^4 - ((T_5 + T_4 + T_3 + T_2)/4)^4)/R_eff_rad;


%% Calculation of Solar absorption by the roof
roof_absorptivity = parameters(6, 9);

% calculate fraction of the day - assumes that t = 0 is midnight
day_num = floor(t/(24*3600)) + 1; % day #
% day 1 is the first day, day 2 is the 2nd day, etc. 

time_frac = t/(24*3600) - (day_num - 1); % fraction of the current day 

q_abs_flux = q_abs(time_frac); % W/m2 

Q_rad_solar_abs = A_roof*q_abs_flux*roof_absorptivity; % [W]

%% Resistances to Heat Transfer

R_5_air = L5/(k5*pi*(r1^2)) + 1/(h1*A_1_bottom);
R_2_air = log(r2/r1)/(2*pi*k2*H_air) + 1/(h1*2*pi*r1*H_air); % K/W
R_air_veg = 2*k_veges/(H_veges*A_1_bottom) + 2*k_air/(H_air*A_1_bottom);
R_2_veg = log(r2/r1)/(2*pi*k2*H_veges) + (2*H2 - H_veges)/(k2*A_1_bottom) + 2*k_veges/(H_veges*A_1_bottom);
R_52 = L5/(k5*pi*(r2^2 - r1^2)) + (H2)/(k2*(pi*(r2^2 - r1^2)));
R_32 = (log(r3/r2)/(2*pi*k3*(H3 - t_s))) + (t_s)/(k3*pi*(r2^2));
R_43 = (log(r4/r3)/(2*pi*k4*dish_thickness)) + (dish_thickness)/(k4*pi*(r3^2));



%% Mass Conservation 


% Lewis number and combined Reynold's analogy 
Le = alpha_inf/D_wa(T_ambient);


% Mass of water evaporated 

MM_H2O = 18.015*10^(-3); % kg/mol


% Capillary Parameters

capillary_vec = parameters(:, 11);
gamma = capillary_vec(1);
d_p = capillary_vec(2);
clay_porosity = capillary_vec(3);
theta_s = capillary_vec(4);
k_clay = capillary_vec(5);

% Evaporative mass flux
m_evap_2 = molar_evap_flux(T_2, T_ambient, RH_inf, h2, rho_inf, Cp_inf, Le)*MM_H2O; % kg/m^2*s

% Effective evaporation area of the clay 
A_c_clay = pi*(r2^2 - r1^2); % m^2 Cross sectional area of clay, perpendicular to z-direction 

% effective Surface area for evaporation based on wetting from capillary
% action
SA_2_eff = 2*pi*r2*h_cap; % [m2]
N_evap_2 = SA_2_eff*m_evap_2/MM_H2O;
N_evap_3 = molar_evap_flux(T_3, T_ambient, RH_inf, h3, rho_inf, Cp_inf, Le)*SA_3;
N_evap_4 = molar_evap_flux(T_4, T_ambient, RH_inf, h4, rho_inf, Cp_inf, Le)*SA_4*0; % mol/s
N_evap_5 = molar_evap_flux(T_5, T_ambient, RH_inf, h5, rho_inf, Cp_inf, Le)*SA_5; % mol/s

total_evap_rate = N_evap_2 + N_evap_3 + N_evap_4 + N_evap_5; % mol/s


%% Evaporative Heat Flows

dH_vap = 2414.3*10^3; % J/kg 

% Heat Flows due to evaporation; All in [W] 
Q_evap_2 = N_evap_2*dH_vap*MM_H2O;
Q_evap_3 = N_evap_3*dH_vap*MM_H2O;
Q_evap_4 = N_evap_4*dH_vap*MM_H2O; 
Q_evap_5 = N_evap_5*dH_vap*MM_H2O;


%% Conductive Heat Flows

Q_5_air = (T_5 - T_air)/R_5_air;
Q_2_air = (T_2 - T_air)/R_2_air; 
Q_air_veges = (T_air - T_veges)/R_air_veg;
Q_2_veg = (T_2 - T_veges)/R_2_veg;
Q_52 = (T_5 - T_2)/R_52;
Q_43 = (T_4 - T_3)/R_43; 
Q_32 = (T_3 - T_2)/R_32; 
  

%% Convective Heat Flow

Q_conv_2 = h2*SA_2*(T_ambient - T_2);
Q_conv_3 = h3*SA_3*(T_ambient - T_3);
Q_conv_4 = h4*SA_4*(T_ambient - T_4);
Q_conv_5 = h5*SA_5*(T_ambient - T_5);
Q_conv_int = h5*A_roof*(T_ambient - T_roof_int);


%% Energy Balance 
% Energy Balance 

dT_airdt = (1/(Cp_air*rho_air*V_air))*(Q_5_air + Q_2_air - Q_air_veges); 
dT_vegdt = (1/(Cp_veges*rho_veges*V_veges))*(Q_air_veges + Q_2_veg);
dT_2dt = (1/(Cp_2*rho_2*V_2))*(Q_52 + Q_32 - Q_2_air + Q_conv_2 + Q_rad_int - Q_evap_2 - Q_2_veg); 
dT_3dt = (1/(Cp_3*rho_3*V_3))*(Q_43 - Q_32 + Q_conv_3 + Q_rad_int - Q_evap_3);
dT_4dt = (1/(Cp_4*rho_4*V_4))*(Q_conv_4 + Q_rad_int - Q_evap_4 ...
    - Q_43); % dish 
dT_5dt = (1/(Cp_5*rho_5*V_5))*(Q_conv_5 + Q_rad_int - Q_evap_5 - Q_5_air - Q_52);
dT_R_ext_dt = (2/(Cp_roof*rho_roof*V_roof)*(-Q_cond_roof + Q_rad_ext + Q_conv_ext + Q_rad_solar_abs));
dT_R_int_dt = (2/(Cp_roof*rho_roof*V_roof)*(Q_cond_roof - Q_rad_int + Q_conv_int));

%% Mass Balance

MM_H2O = 18.015; % g/mol 
rho_water = 1000; % kg/m3 


if V_water <= 0
    dV_dt = 0; 
    % infinite water in device; otherwise, set evaporative flows to 0 
else
    dV_dt = -(total_evap_rate*MM_H2O)/(1000*rho_water); % m3/s 
end
 

%% Capillary Effects - Mass and Momentum balance  

% Properties and Constants
mu_H2O = 8.9E-4; % Pa*s; water viscosity at 25C 
rho_H2O = 1000; % kg/m3 
g = 9.8; % [m/s^2]

% Capillary Pressure

P_c = gamma*cos(theta_s)/(d_p/2); % Pa 

% Gravitational Pressure 
P_g = rho_H2O*g*h_cap; 

% Viscous Pressure Terms 

% Due to evaporation

P_v = (h_cap^2)*r2*pi*m_evap_2*mu_H2O/(k_clay*rho_H2O*A_c_clay);


% Gov equation based on mass and momentum balance 

dhdt = (P_c - P_g - P_v)*k_clay/(clay_porosity*mu_H2O*h_cap); % [m/s] 


%% Derivatives
derivs(1) = dT_airdt;
derivs(2) = dT_vegdt;
derivs(3) = dT_2dt;
derivs(4) = dT_3dt;
derivs(5) = dT_4dt;
derivs(6) = dT_5dt;
derivs(7) = dV_dt;
derivs(8) = dT_R_ext_dt;
derivs(9) = dT_R_int_dt;
derivs(10) = dhdt;


end