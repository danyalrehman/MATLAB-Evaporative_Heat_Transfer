function [derivs] = Empty_Pot_ODEs(t, vars, parameters, other)


derivs = zeros(length(vars), 1);

%% Define Variables 

% All in Kelvin 


% Temperature 
T_1 = vars(1);% inner chamber, just air 
T_2 = vars(2); 
T_3 = vars(3);
T_4 = vars(4);
T_5 = vars(5);
V_water = vars(6); % m3 of water in device 

T_roof_ext = vars(7); % Kelvin, external roof surface temp
T_roof_int = vars(8); % Kelvin, intenral roof surface temp 


% unpack the ambient data
t_data = other(:, 1); % seconds; time data 
T_outside = other(:, 2) + 273.15; % Kelvin; Temperature Data for ambient 
outside_humidity = other(:, 3); % Ambient humidity data 
air_props = other(:, 4:7); % Calculated air properties 

mu_inf = interp1(t_data, air_props(:,1), t); % viscosity of ambient air 
k_inf = interp1(t_data, air_props(:,2), t); % thermal conductivity of ambient air 
Cp_inf = interp1(t_data, air_props(:,3), t); % Cp of ambient air 
rho_inf = interp1(t_data, air_props(:,4), t); % density of ambient air 

nu_inf = mu_inf/rho_inf; % kinematic viscosity of ambient air 

alpha_inf = k_inf/(rho_inf*Cp_inf); % thermal diffusivity of ambient air 

T_ambient = interp1(t_data, T_outside, t); % T_inf for time step t 

RH_inf = interp1(t_data, outside_humidity, t); % Relative humidity at time step t 


%% Unpack all of the parameters 
% Volumes
V_vec = parameters(:,3); % Volume 

V_1 = V_vec(1);
V_2 = V_vec(3);
V_3 = V_vec(4);
V_4 = V_vec(5);
V_5 = V_vec(6);

% Unpack all the parameters 
Cp_vec = parameters(:, 2); % Cp

Cp_1 = Cp_sat_air(T_1);
Cp_2 = Cp_vec(3);
Cp_3 = Cp_vec(4); % Sand 
Cp_4 = Cp_vec(5);
Cp_5 = Cp_vec(6);

Cp_water = 4184; % J/kgK
x_water = V_water/V_3; 
C_H2O_3 = (V_water*1E6)/(V_3*18.015); % moles of water/m3 sand 
C_H2O_inf = C_water_IG(T_ambient, RH_inf);  % mol/m3 


% sand_porosity = parameters(5, 9);

% Weighted average 
Cp_3 = Cp_water*x_water + (1-x_water)*Cp_3;

outer_speed = parameters(2, 9); 
inner_speed = parameters(1, 9); 
speeds = [inner_speed, outer_speed];


% Densities 
rho_vec = parameters(:,1); % Density

rho_1 = rho_sat_air(T_1); 
rho_2 = rho_vec(3);
rho_3 = rho_vec(4);
rho_4 = rho_vec(5);
rho_5 = rho_vec(6);

% Thermal Conductivity 
k_vec = parameters(:, 6); % W/mK
k1 = k_sat_air(T_1);
k2 = k_vec(3); % inner pot 
k3 = k_vec(4);  % sand
k4 = k_vec(5); % outer pot 
k5 = k_vec(6); % Dish

% Radii and Heights 

r_vec = parameters(:, 4);
r1 = r_vec(2); 
r2 = r_vec(3);
r3 = r_vec(4);
r4 = r_vec(5);

H_vec = parameters(:, 5);

H1 = H_vec(1);
H2 = H_vec(3);
H3 = H_vec(4);
H4 = H_vec(5);
L5 = H_vec(6);

% Heat Transfer Coefficients 
ambient_vars = [T_ambient, RH_inf, mu_inf, k_inf, Cp_inf, rho_inf];

h_vec = ConvectionModel(speeds, vars, t, parameters, ambient_vars);
h1 = h_vec(1);
h4 = h_vec(2);
h5 = h_vec(3);

% Cloud cover factor  (Clark and Allen)
N = 0; % Clear skies 
CC = 1 + 0.0224*N - 0.0035*N^2 + 0.00028*N^3;  


% Calculate T_sky for long-wave radiation
% Calcualte Ambeint Dew Point Temperature
Tdp_amb = T_dp(T_ambient, RH_inf);
Tsky = T_sky(T_ambient, E_clear(Tdp_amb), CC);

E_clay = parameters(1, 8);
E_cloth = parameters(2, 8); 
sigma = 5.6703*10^(-8); % W/m^2K^4


%% Short wave radiation - Solar irradiance 
roof_absorptivity = parameters(5, 8);

% calculate what fraction of the day it is - assumes that t = 0 is midnight
day_num = floor(t/(24*3600)) + 1; % day #
% day 1 is the first day, day 2 is the 2nd day, etc. 

time_frac = t/(24*3600) - (day_num - 1); % fraction of the current day 

q_abs_flux = q_abs(time_frac); % W/m2; Incident radiation 

%%  Surface Areas
SA_vec = parameters(:,7);
A_1_bottom = SA_vec(2); % Area of contact between bottom of surface 1 and contacting part of
                        % surface 2 
SA_4 = SA_vec(5);
SA_5 = SA_vec(6);

%% Radiation Networks and Roof calculations 


z_roof = parameters(1, 10); % [m]; Height from ground to the roof 
house_width = z_roof; % [m]
house_length = house_width; % [m] 

u_roof = parameters(2, 10); % [m/s]
A_roof = house_width*house_length; % m2

Re_roof_ext = rho_inf*u_roof*(sqrt(A_roof))/(mu_inf);
Pr_inf = nu_inf/alpha_inf;
% Assume inner ceiling wind speed to be 0.5 m/s
Re_roof_int = rho_inf*0.5*(sqrt(A_roof))/(mu_inf);


A_wall = z_roof*house_width; % [m2]
A_g = house_width*house_length - SA_5; % Area of the ground [m2]

E_wall = parameters(4,8);

% Assume only forced 
Nu_ext_roof = 0.664*sqrt(Re_roof_ext)*(Pr_inf^(1/3));
Nu_int_roof = 0.664*sqrt(Re_roof_int)*(Pr_inf^(1/3));

% Ra_roof = g*B_outside*(T_roof_ext - Tsky)*(L_roof^3)/(nu_inf*alpha_inf);
h_conv_roof_ext = Nu_ext_roof*(k_inf)/(sqrt(A_roof)); % [W/m2K]
h_conv_roof_int = Nu_int_roof*(k_inf)/(sqrt(A_roof)); % [W/m2K]

% Forced Convection above roof 
E_roof = parameters(3,8); 


h_rad_roof_ext = E_roof*sigma*(T_roof_ext + Tsky)*(T_roof_ext^2 + Tsky^2);  

% Conduction in the roof 
L_roof = parameters(3, 10); % [m]  thickness of the roof



k_roof = parameters(6, 10); % [W/mK] 
rho_roof = parameters(4, 10); % kg/m3
Cp_roof = parameters(5, 10); % [J/kgK]

Q_cond_roof = k_roof*A_roof*(T_roof_ext - T_roof_int)/L_roof;
Q_rad_ext = h_rad_roof_ext*A_roof*(Tsky - T_roof_ext); % long wave, [W]
Q_rad_solar_abs = A_roof*q_abs_flux*roof_absorptivity; % [W]
Q_conv_ext = h_conv_roof_ext*A_roof*(T_ambient - T_roof_ext);

V_roof = L_roof*A_roof; % [m3] 

% Interior

% Surface 0 represents the inside walls and roof of the dwelling 

E_0 = (E_roof + E_wall)/2;
E_s = (E_clay + E_cloth)/2; 

A_0 = A_roof + 4*A_wall; 
A_s = SA_4 + SA_5; 

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

Q_rad_int = sigma*(((T_roof_int + T_ambient)/2)^4 - ((T_5 + T_4)/2)^4)/R_eff_rad; % [W] 


%% Resistances 
% Heat Transfer
 
R_51 = L5/(k5*pi*(r1^2)) + 1/(h1*pi*(r1^2));
R_21 = log(r2/r1)/(2*pi*k2*2*H2) + (H2-H1)/(k2*A_1_bottom) + 1/(h1*A_1_bottom); % K/W
R_52 = L5/(k5*pi*(r2^2 - r1^2)) + (H2)/(k2*(pi*(r2^2 - r1^2)));
R_32 = log(r3/r2)/(2*pi*k3*2*H2) + (H3-H2)/(k3*pi*(r2^2));
R_53 = L5/(k5*pi*(r3^2 - r2^2)) + (H3)/(k3*(pi*(r3^2 - r2^2)));
R_43 = log(r4/r3)/(2*pi*k4*2*H3) + (H4-H3)/(k4*pi*(r3^2));
R_54 = L5/(k5*pi*(r4^2 - r3^2)) + (H4)/(k4*pi*(r4^2 - r3^2));


%% Mass Conservation  

% Lewis number and combined Reynold's analogy 
Le = alpha_inf/D_wa(T_ambient); 

h_m_4 = h4/(rho_inf*Cp_inf*Le^(2/3));
h_m_5 = h5/(rho_inf*Cp_inf*Le^(2/3));


C_H2O_4_max = C_water_IG(T_4, 1); % thin saturated film 
C_H2O_5_max = C_water_IG(T_5, 1); % thin saturated film 

N_evap_4 = molar_evap_flux(C_H2O_4_max, C_H2O_inf, h_m_4)*SA_4; % mol/s 

N_evap_5 = molar_evap_flux(C_H2O_5_max, C_H2O_inf, h_m_5)*SA_5; % mol/s

total_evap_rate = N_evap_4 + N_evap_5; % mol/s
 
%% Evaporative Heat Flows

dH_vap = 2414.3*10^3; % J/kg 

MM_H2O = 18.015*10^(-3); % kg/mol 

Q_evap_4 = N_evap_4*MM_H2O*dH_vap; % [W]
Q_evap_5 = N_evap_5*MM_H2O*dH_vap; % [W]


%% Conductive Heat Flows [Watts] 

Q_51 = (T_5 - T_1)/R_51;
Q_21 = (T_2 - T_1)/R_21; 
Q_52 = (T_5 - T_2)/R_52;
Q_53 = (T_3 - T_2)/R_53; 
Q_43 = (T_4 - T_3)/R_43; 
Q_32 = (T_3 - T_2)/R_32; 
Q_54 = (T_5 - T_4)/R_54;


  
%% Convective Heat Flow [Watts]
Q_conv_4 = h4*SA_4*(T_ambient - T_4); 
Q_conv_5 = h5*SA_5*(T_ambient - T_5);
Q_conv_int = h_conv_roof_int*A_roof*(T_ambient - T_roof_int);


%% Mass Balance

MM_H2O = 18.015; % g/mol 
rho_water = 1000; % kg/m3 

if V_water <= 0
    dV_dt = 0; 
  
else
    dV_dt = -(total_evap_rate*MM_H2O)/(1000*rho_water); % m3/s 
end



%% Energy Balance 
% Energy Balance 

dT_1dt = (1/(Cp_1*rho_1*V_1))*(Q_51 + Q_21); 
dT_2dt = (1/(Cp_2*rho_2*V_2))*(Q_52 + Q_32 - Q_21); 
dT_3dt = (1/(Cp_3*rho_3*V_3))*(Q_53 + Q_43 - Q_32);
dT_4dt = (1/(Cp_4*rho_4*V_4))*(Q_54 + Q_conv_4 + Q_rad_int - Q_evap_4 ...
    - Q_43);
dT_5dt = (1/(Cp_5*rho_5*V_5))*(Q_conv_5 + Q_rad_int - Q_evap_5 - Q_51 ...
    - Q_52 - Q_53 - Q_54);
dT_R_ext_dt = (2/(Cp_roof*rho_roof*V_roof)*(-Q_cond_roof + Q_rad_ext + Q_conv_ext + Q_rad_solar_abs));
dT_R_int_dt = (2/(Cp_roof*rho_roof*V_roof)*(Q_cond_roof - Q_rad_int + Q_conv_int));


%% Derivatives

derivs(1) = dT_1dt;
derivs(2) = dT_2dt;
derivs(3) = dT_3dt;
derivs(4) = dT_4dt;
derivs(5) = dT_5dt;
derivs(6) = dV_dt; 
derivs(7) = dT_R_ext_dt;
derivs(8) = dT_R_int_dt; 


end