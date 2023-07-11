
% Script for Pot in Dish modeling 
% Authors: Ethan McGarrigle, Danyal Rehman; Massachusetts Institute of
% Technology 

% "Run" this script 

clear all, close all, clc
tic
%% Load any experimental data (Ambient inputs, observed device temperature)
 
data = readmatrix('Sample_Experimental_Data.xls');

t_recorded = data(1:end, 1);
T_inside = data(1:end, 2);
T_outside = data(1:end, 3);
outside_humidity = data(1:end, 4);

% Calculate the air fluid properties as a function of ambient temperature and
% humidity 
air_props = calcAirProps(T_outside, outside_humidity);


%% Ambient Conditions and Inputs
RH_inf_initial = outside_humidity(1); 
T_1_initial = T_inside(1) + 273.15; % Kelvin

t_end = t_recorded(end); % [s]

% Specify air speeds 
outer_speed = 0.5; % velocity in m/s;  
inner_speed = 0;   % velocity in m/s;  

%% Geometric inputs and Material inputs

% All units are SI 


% Vegetables: mass and properties


mass_tomato = 2; % kg
rho_tomato = 1022.87;  % kg/m3

V_veges = mass_tomato/rho_tomato;
k_veges = 0.59; % W/mK (Tomato) 
rho_veges = rho_tomato;


% Geometric Inputs

sand_thickness = 5*10^(-2); % [m] Thickness of the sand 
pot_thickness = 2*10^(-2); % [m] Thickness of the clay pot 
dish_thickness = 3E-3; % [m] Thickness of the plastic dish 

r1 = 0.1670; % [m] 
r2 = r1 + pot_thickness;
r3 = r2 + sand_thickness;
r4 = r3 + dish_thickness; 


% Heights 

H4 = 0.2; % [m] Dish height
H3 = H4 - dish_thickness; % Sand height 
H2 = 0.38; % height of the full pot 
H1 = H2 - pot_thickness; 


H_veges = V_veges/(pi*r1^2); 
H_air = H1 - H_veges; 

t_s = (3.5E-2); % [m], specified sand depth 


L5 = 0.0025; % thickness of lid/cloth, meters

H_i = [H_air; H_veges; H2; H3; H4; L5];

V_1 = H1*pi*(r1^2); % m3 % Volume of air in inner chamber 
V_air = V_1 - V_veges; 
V_2 = pi*H2*r2^2 - V_1; % Volume of inner pot
V_3 = pi*(r3^2)*t_s + pi*(r3^2 - r2^2)*(H3 - t_s);
V_4 = pi*(r4^2)*(dish_thickness) + pi*(r4^2 - r3^2)*H3; 
V_5 = L5*pi*r2^2; % Volume of lid/cloth 

V_i = [V_air; V_veges; V_2; V_3; V_4; V_5];

if V_air < 0
   warning('Inner chamber not large enough to fill specified vegetable volume. Results may be inaccurate');   
end

k_air = k_sat_air(T_1_initial); % W/mK
k1 = 0.0265;  % air at 100% humidity, 30C  
k2 = 1.0; % [W/mK] (clay)
k3 = 3.27; % [W/mK] (sand)         
k4 = 0.16; % [W/mK] %  dish
k5 = 3; % Cloth or lid thermal conductivity 

k_i = [k_air; k_veges; k2; k3; k4; k5];

sand_porosity = 0.40; 

r_i = [NaN; r1; r2; r3; r4; t_s];


% Tomatoes are 94% Water and 6% cellulose 
Cp_water = 4184; % J/kgK
Cp_cellulose = 1400; % J/kgK 
Cp_veg = (1720 + 960)/2; % J/kgK
Cp_tomato = 0.06*(Cp_cellulose) + Cp_water*0.94; % J/kgK

Cp_2 = 900; % J/kgK % Pot
Cp_3 = 1532.7; % 
Cp_4 = 2800; % J/kgK % dish  
Cp_5 = 1360; % J/kgK % Cloth? 

Cp_i = [Cp_sat_air(T_1_initial); Cp_veg; Cp_2; Cp_3; Cp_4; Cp_5];

rho_air_sat = 1.18; % kg/m3
rho_tomato = 1022.87;  % kg/m3
rho_1 = rho_air_sat;  % kg/m3
rho_2 = 1900; % kg/m3
rho_3 = 2057; % kg/m3 
rho_4 = 946; % kg/m3 
rho_5 = 1.46*1000; % kg/m3; Cloth density 

rho_i = [rho_sat_air(T_1_initial); rho_veges; rho_2; rho_3; rho_4; rho_5]; 

% Area of contact between bottom of surface 1 and contacting part of
% surface 2 
A_1_bottom = pi*(r1^2); 

SA_2 = 2*pi*(H2 - (H3 - t_s))*r2;
SA_3 = pi*(r3^2 - r2^2);
SA_4 = pi*r4*2*H4; % m2
SA_5 = L5*2*pi*r2 + pi*r2^2; % m2

SA_i = [NaN; A_1_bottom; SA_2; SA_3; SA_4; SA_5];

% Emissivity
E_clay = 0.775; % Emissivity of clay 
E_cloth = 0.77; % Emissiviy of cloth
E_dish = 0.95; % Dish emissivity 
E_sand = 0.76; % Sand 

V_capacity = V_3*sand_porosity; % m3


% Put all of the properties in a vector 
geo_and_props = [rho_i, Cp_i, V_i, r_i, H_i, k_i, SA_i];


E_wall = 0.7; % Emissivity of the wall (clay?) 


z_roof = 5; % [meters]
L_roof = 0.02;

% Assume an aluminum roof 
k_roof = 168; % [W/mK] 
rho_roof = 2790; % kg/m3
Cp_roof = 883; % [J/kgK]
roof_solar_absorptivity = 0.08; 
E_roof = 0.02; 

u_roof = 5; % [m/s]

roof_vec = [z_roof; u_roof; L_roof; rho_roof; Cp_roof; k_roof];


% Capillary Effects
gamma = 0.0728; % [N/m] Surface tension between water and air 
d_p = 0.6385E-6; % [m] effective pore diameter of clay; avg from the table
clay_porosity = 0.3348; % average value 
theta_s = 0; % Contact Angle in radians 
k_clay = (d_p^2)*(clay_porosity^3)/(150*(1 - clay_porosity)^2); % m^2 
capillary_vec = [gamma; d_p; clay_porosity; theta_s; k_clay; NaN];


 

%% Set up solver (initial conditions, time domain, etc.)


capillary_height_i = 0.60*(H2 - H4); % [m] % Initial height of capillary rise - approximate Steady State value from observation
T_roof_ext_i = T_outside(1) + 273.15;
T_roof_int_i = T_outside(1) + 273.15;

% If the pot is empty, use the appropriate function 
if V_veges == 0
    f_handle_Transient = @ODE_fxn_Empty_Pot; % empty inner chamber
    initial_cond = [T_1_initial, T_1_initial, T_1_initial, T_1_initial , ...
       T_1_initial,  V_capacity, T_roof_ext_i, T_roof_int_i, capillary_height_i]; 
else
    f_handle_Transient = @ODE_fxn_Filled_Pot; % Vegetables are present
    initial_cond = [T_1_initial, T_1_initial, T_1_initial, T_1_initial, T_1_initial , ...
       T_1_initial,  V_capacity, T_roof_ext_i, T_roof_int_i, capillary_height_i]; 
end


constraints_vec = [inner_speed; outer_speed; NaN; V_capacity; sand_porosity; roof_solar_absorptivity];
emissivity_vec = [E_sand; E_dish; E_clay; E_cloth; E_roof; E_wall]; 
parameters = [geo_and_props, emissivity_vec, constraints_vec, roof_vec, capillary_vec];

t_span = [0, t_recorded(end)];

other = [t_recorded, T_outside, outside_humidity, air_props, T_inside];

%% Solve 

% Solve using ode23s; t is the time in seconds; sol represents the solution
% vector
[t, sol] = ode23s(f_handle_Transient, t_span, initial_cond, [], parameters, other);


%% Post-model Processing and Plotting

% Recalcualte the wet bulb temperatures
Twb_vecinterp = zeros(length(t), 1);
Twb_obs = zeros(length(t_recorded), 1);

for t_step = 1:length(t)
    time = t(t_step); % seconds 
    T_ambient = interp1(t_recorded, T_outside, time);
    humidity_ambient = interp1(t_recorded, outside_humidity, time); 
    [Tdb, w, phi, h, Tdp, v, Twb_vecinterp(t_step)] = Psychrometricsnew('Tdb', T_ambient, ...
         'phi', humidity_ambient*100);
end

for t_step = 1:length(t_recorded)
    time = t_recorded(t_step); % seconds 
    T_ambient = T_outside(t_step);
    humidity_ambient = outside_humidity(t_step); 
    [Tdb, w, phi, h, Tdp, v, Twb_obs(t_step)] = Psychrometricsnew('Tdb', T_ambient, ...
         'phi', humidity_ambient*100);
end

H2O_flux_2 = [];
H2O_flux_3 = [];
H2O_flux_4 = [];
H2O_flux_5 = [];
Le = [];
molar_evap_rate_2 = [];
%  Capillary rise height is column 9 if no vegetable mass, col 10 otherwise 
if V_veges == 0
    SA_2_eff = 2.*pi.*r2.*sol(:,9); % [m2] 
else
    SA_2_eff = 2.*pi.*r2.*sol(:,10); % [m2] 
end


for i = 1:length(t)

    T_ambient = interp1(t_recorded, T_outside, t(i)) + 273.15; % Kelvin
    humidity_ambient = interp1(t_recorded, outside_humidity, t(i)); 
    rho_inf = interp1(t_recorded, air_props(:,4), t(i));
    Cp_inf = interp1(t_recorded, air_props(:,3), t(i));
    mu_inf = interp1(t_recorded, air_props(:,1), t(i));
    k_inf = interp1(t_recorded, air_props(:,2), t(i));
    alpha_inf = k_inf/(rho_inf*Cp_inf); % m2/s 
    
    Le(i) = alpha_inf/D_wa(T_ambient);
    speeds_vec = [inner_speed, outer_speed];
    
    ambient_vars = [T_ambient, humidity_ambient, mu_inf, k_inf, Cp_inf, rho_inf];
    h_vec = ConvectionModel(speeds_vec, sol(i, :), t(i), parameters, ambient_vars);
    if V_veges == 0
    H2O_flux_2(i) = molar_evap_flux(sol(i,2), T_ambient, humidity_ambient, ...
        h_vec(2), rho_inf, Cp_inf, Le(i)); % mol/m2*s
    H2O_flux_3(i) = molar_evap_flux(sol(i,3), T_ambient, humidity_ambient, ...
        h_vec(3), rho_inf, Cp_inf, Le(i)); % mol/m2*s
    H2O_flux_4(i) = molar_evap_flux(sol(i,4), T_ambient, humidity_ambient, ...
        h_vec(4), rho_inf, Cp_inf, Le(i))*0; % mol/m2*s % no evaporation off of dish 
    H2O_flux_5(i) = molar_evap_flux(sol(i,5), T_ambient, humidity_ambient, ...
        h_vec(5), rho_inf, Cp_inf, Le(i)); % mol/m2*s
    else
        H2O_flux_2(i) = molar_evap_flux(sol(i,3), T_ambient, humidity_ambient, ...
            h_vec(2), rho_inf, Cp_inf, Le(i)); % mol/m2*s
        H2O_flux_3(i) = molar_evap_flux(sol(i,4), T_ambient, humidity_ambient, ...
            h_vec(3), rho_inf, Cp_inf, Le(i)); % mol/m2*s
        H2O_flux_4(i) = molar_evap_flux(sol(i,5), T_ambient, humidity_ambient, ...
            h_vec(4), rho_inf, Cp_inf, Le(i))*0; % mol/m2*s % no evaporation off of dish 
        H2O_flux_5(i) = molar_evap_flux(sol(i,6), T_ambient, humidity_ambient, ...
            h_vec(5), rho_inf, Cp_inf, Le(i)); % mol/m2*s
    end
    
    
    molar_evap_rate_2(i) = H2O_flux_2(i)*SA_2_eff(i);  % mol/s
    
end

molar_H2O_evap_rate = molar_evap_rate_2 + H2O_flux_3*SA_3 + H2O_flux_4.*SA_4.*0 + H2O_flux_5.*SA_5; % mol/s

MM_H2O = 18.015*10^(-3); % kg/mol 
rho_water = 1000; % kg/m3

volume_evap_rate = molar_H2O_evap_rate.*MM_H2O./rho_water; % m3/s 
% Convert from m3 to mL 
volume_evap_rate = volume_evap_rate.*1000; % L/s 

% Integrate using Trapz to get the total volume lost 
total_water_lost = trapz(t, volume_evap_rate); % L 

disp(['Total Water Evaporated: ', num2str(total_water_lost), ' [L]']); 

Twb_vecinterp = Twb_vecinterp(1:length(Twb_vecinterp));

% Refrigeration Efficiency 
eta_model = ((sol(:,1) - 273.15) - interp1(t_recorded, T_outside, t))./(Twb_vecinterp - interp1(t_recorded, T_outside, t)); 
eta_obs = (T_inside - T_outside)./(Twb_obs - T_outside); 


figure(1); 
hold on;
plot(t./3600, sol(:,1) - 273.15, '-', 'linewidth', 1.5);
% % % Plot ambient and wetbulb temp
plot(t./3600, interp1(t_recorded, T_outside, t), '-', 'linewidth', 2.5); % ambient 
plot(t_recorded(1:(end))./3600, T_inside(1:(end)), 'linewidth', 2.5); % observed 
plot(t./3600, Twb_vecinterp, '-', 'linewidth', 1.5);
legend('Modeled', 'T_{\infty}', 'Observed', 'T_{wb}');
xlabel('\textbf{Time [hours]}', 'interpreter', 'latex');
ylabel('\textbf{Temperature [C]}', 'interpreter', 'latex');  
set(gca,'FontSize',16,'LineWidth',2,'FontName','Arial');
set(gcf,'Color','w');
figure(1);

T_interp_model = interp1(t, sol(:,1) - 273.15, t_recorded); % modeled data


%% Error Estimation
A = isnan(T_interp_model);
nan_index_array = []; 
k = 1; 

for i = 1:length(T_interp_model)
   if A(i) == 1
      nan_index_array(k) = i; 
      k = k + 1;    
   end 
end

% loop through all of T_inside and checks the indices 
% shorten T_inside via the loop 
for i = 1:length(T_inside)  
   for k = 1:length(nan_index_array) 
       if i == nan_index_array(k)
           T_inside(i) = NaN; 
       end
       
   end
    
end


T_interp_model = rmmissing(T_interp_model); % modeled data 
T_inside = rmmissing(T_inside); % obseved 


% Error Calculations
T_obs_mean = mean(T_inside); % experimental 
SS_total = sum((T_inside - T_obs_mean).^2); % experimental compared to experimental mean 
SS_res = sum((T_inside -  T_interp_model).^2);

T_model_mean = mean(sol(:, 1) - 273.15);

NMSE = (SS_res/(T_model_mean*T_obs_mean))/(length(sol(:,1))); % normalized mean squared error

MAE = sum(abs(T_interp_model - T_inside))/length(T_interp_model); % Mean absolute error 

RMSE = (SS_res/(length(T_inside)))^(1/2); % root mean squared error 

toc





