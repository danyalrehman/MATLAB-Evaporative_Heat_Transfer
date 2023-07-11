
% Script for Pot in Pot modeling 
% Authors: Ethan McGarrigle, Danyal Rehman; Massachusetts Institute of
% Technology 

% "Run" this script 

clear all, close all, clc
tic

%% Read/Load the experimental data 

data = xlsread('Sample_Experimental_Data.xls');

t_recorded = data(1:end, 1);
T_inside = data(1:end, 2);
T_outside = data(1:end, 3);
outside_humidity = data(1:end, 4);

% Calculate the air fluid properties as a function of ambient temperature and
% humidity 
air_props = calcAirProps(T_outside, outside_humidity);


%% Ambient Conditions and Environmental Inputs
RH_inf_initial = outside_humidity(1); 

% Specify air speeds 
outer_speed = 0.5; % velocity in m/s;  
inner_speed = 0;   % velocity in m/s;  

%% Geometric inputs and Material inputs

% All units are SI 


% Vegetables: mass and properties

mass_tomato = 0; % kg
rho_tomato = 1022.87;  % kg/m3

V_veges = mass_tomato/rho_tomato;
k_veges = 0.59; % W/mK (Tomato) 
rho_veges = rho_tomato;

sand_thickness = 4*10^(-2); % [m] Thickness of the sand  
pot_thickness = 1.5*10^(-2); % Thickness of the pots


% Region 1 is the inner chamber consisting of vegetables and humid air
% Region 2 is Inner Pot
% Region 3 is sand 
% Region 4 is Outer Pot
% Region 5 is the cloth on top of the device 


% CYLINDRICAL COORDINATES 

% Geometry

% r = 0 occurs at the center of the device, in the inner chamber 
 
r1 = 0.1050; % [m] (base case) 
r2 = r1 + pot_thickness; 
r3 = r2 + sand_thickness;
r4 = r3 + pot_thickness;


r_i = [NaN; r1; r2; r3; r4; NaN];

% Heights; Defined as Half-heights

H4 = (0.3/2); % [m] half height of outer clay pot 
H3 = H4 - pot_thickness; % half height of sand layer 
H2 = H3 - sand_thickness; % half height of inner pot 
H1 = H2 - pot_thickness; % half height of inner chamber 

V_1 = 2*H1*pi*(r1^2); % m3 % Total volume of inner chamber

% Vegetable and Air heights are NOT defined as half-heights
H_veges = V_veges/(2*pi*r1^2); % half height of vegetable mass in inner chamber 
H_air = (2*H1 - 2*H_veges)/2; % half height of air in inner chamber 

L5 = 0.0025; % [m] thickness of lid/cloth 

H_i = [H_air; H_veges; H2; H3; H4; L5];


V_air = V_1 - V_veges; % Volume of air [m3]
V_2 = 2*pi*H2*r2^2 - V_1; % Volume of inner pot
V_3 = 2*pi*H3*r3^2 - (V_2 + V_1); % Volume of Sand
V_4 = 2*pi*H4*r4^2 - (V_2 + V_1 + V_3); % Volume of outer clay 
V_5 = L5*pi*r4^2; % Volume of lid/cloth 

V_i = [V_air; V_veges; V_2; V_3; V_4; V_5];

if V_air < 0
   warning('Inner chamber not large enough to fill specified vegetable volume. Results may be inaccurate');   
end

T_1_initial = T_inside(1) + 273.15;

k1 = k_sat_air(T_1_initial); % thermal conductivity of air at 100% Humidity  
k2 = 1.3; % [W/mK] 
k3 = 3.27; % [W/mK] 
k4 = 2.0; % [W/mK] 
k5 = 0.44; %  Cloth or lid thermal conductivity 

k_air = k_sat_air(T_1_initial); % W/mK

k_i = [k_air; k_veges; k2; k3; k4; k5];

sand_porosity = 0.40; 

Cp_air = Cp_sat_air(T_1_initial); % J/kgK
% Assumes tomatoes are 94% Water and 6% cellulose 
Cp_water = 4184; % J/kgK
Cp_cellulose = 1400; % J/kgK 
Cp_tomato = 0.06*(Cp_cellulose) + Cp_water*0.94; % J/kgK

Cp_veg = Cp_tomato; % J/kgK 

Cp_1 = Cp_air; % J/kgK 
Cp_2 = 900; % J /kgK
Cp_3 = 1532.7; % J/kgK
Cp_4 = 2423; % J/kgK % Pot 
Cp_5 = 1360; % J/kgK % Cloth

Cp_i = [Cp_sat_air(T_1_initial); Cp_veg; Cp_2; Cp_3; Cp_4; Cp_5];


rho_2 = 2250; % kg/m3
rho_3 = 2057; % kg/m3 
rho_4 = 2250; % kg/m3
rho_5 = 1.46*1000; % kg/m3

rho_i = [rho_sat_air(T_1_initial); rho_veges; rho_2; rho_3; rho_4; rho_5]; 


% Area of contact between bottom of surface 1 and contacting part of
% surface 2 (Surface Area of the bottom of the inner chamber) 
A_1_bottom = pi*(r1^2); % m2 

SA_4 = 2*pi*r4*2*H4; % m2; Outer surface area of the clay pot, in contact with the air 
SA_5 = L5*2*pi*r4 + pi*r4^2; % m2; Outer surface area of the cloth/lid, in contact with the air 

SA_i = [NaN; A_1_bottom; NaN; NaN; SA_4; SA_5];

E_clay = 0.775; % Emissivity of clay 
E_cloth = 0.77; % Emissiviy of cloth; 

V_H2O_i = V_3*sand_porosity; % m3; Initial Volume of water in saturated sand 


% Put all of the properties in a vector 
geo_and_props = [rho_i, Cp_i, V_i, r_i, H_i, k_i, SA_i];

E_wall = 0.7; % Emissivity of the dwelling walls (assumed clay) 


z_roof = 5; % [meters] 
L_roof = 0.02; % [m]

% Assume an aluminum roof 
k_roof = 168; % [W/mK] aluminum
rho_roof = 2790; % [kg/m3] aluminum
Cp_roof = 883; % [J/kgK] aluminum
roof_solar_absorptivity = 0.08; 
E_roof = 0.02; 

u_roof = 5; % [m/s]

roof_vec = [z_roof; u_roof; L_roof; rho_roof; Cp_roof; k_roof];


%% Set Initial Conditions, initial solver, set up time interval to run the simulation 

t_span = [0, t_recorded(end)];

% Set initial conditions 
T_roof_ext_i = T_outside(1) + 273.15;
T_roof_int_i = T_outside(1) + 273.15;

if V_veges == 0
    f_handle_Transient = @Empty_Pot_ODEs; % empty inner chamber
    initial_cond = [T_1_initial, T_1_initial, T_1_initial, T_1_initial , ...
    T_1_initial,  V_H2O_i, T_roof_ext_i, T_roof_int_i]; 
else
    f_handle_Transient = @Filled_Pot_ODEs; 
    initial_cond = [T_1_initial, T_1_initial, T_1_initial, T_1_initial, T_1_initial , ...
    T_1_initial,  V_H2O_i, T_roof_ext_i, T_roof_int_i]; 
end



constraints_vec = [inner_speed; outer_speed; NaN; V_H2O_i; sand_porosity; NaN];
emissivity_vec = [E_clay; E_cloth; E_roof; E_wall; roof_solar_absorptivity; NaN];
parameters = [geo_and_props, emissivity_vec, constraints_vec, roof_vec];

% Put in the recorded data and air propeties for the experiment
other = [t_recorded, T_outside, outside_humidity, air_props, T_inside];

%% Run the simulation

% Solve using ode23s; t is the time in seconds; sol represents the solution
% vector
[t, sol] = ode23(f_handle_Transient, t_span, initial_cond, [], parameters, other);


%% Wet Bulb

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

 
H2O_flux_4 = [];
H2O_flux_5 = [];
Le = [];
h_m_4 = [];
h_m_5 = [];
T_drop = [];

% Evaluate the model at the solution 
for i = 1:length(t)
    % Calculate convection coefficients 
    T_ambient = interp1(t_recorded, T_outside, t(i)) + 273.15; % Kelvin
    humidity_ambient = interp1(t_recorded, outside_humidity, t(i)); 
    rho_inf = interp1(t_recorded, air_props(:,4), t(i));
    Cp_inf = interp1(t_recorded, air_props(:,3), t(i));
    mu_inf = interp1(t_recorded, air_props(:,1), t(i));
    k_inf = interp1(t_recorded, air_props(:,2), t(i));
    % Lewis number and combined Reynold's analogy 
    alpha_inf = k_inf/(rho_inf*Cp_inf); % m2/s 
    
    Le(i) = alpha_inf/D_wa(T_ambient);
    
    speeds = [inner_speed, outer_speed];
    ambient_vars = [T_ambient, humidity_ambient, mu_inf, k_inf, Cp_inf, rho_inf];
    h_vec = ConvectionModel(speeds, sol(i, :), t(i), parameters, ambient_vars);
    h_m_4(i) = h_vec(2)./(rho_inf.*Cp_inf*(Le(i).^(2/3)));
    h_m_5(i) = h_vec(3)./(rho_inf.*Cp_inf.*(Le(i).^(2/3)));
    if V_veges == 0
       H2O_flux_4(i) = molar_evap_flux(C_water_IG(sol(i, 4), 1), C_water_IG(T_ambient, humidity_ambient), h_m_4(i)); % mol/m2*s
       H2O_flux_5(i) = molar_evap_flux(C_water_IG(sol(i, 5), 1), C_water_IG(T_ambient, humidity_ambient), h_m_5(i)); % mol/m2*s 
    else
       H2O_flux_4(i) = molar_evap_flux(C_water_IG(sol(i, 5), 1), C_water_IG(T_ambient, humidity_ambient), h_m_4(i)); % mol/m2*s
       H2O_flux_5(i) = molar_evap_flux(C_water_IG(sol(i, 6), 1), C_water_IG(T_ambient, humidity_ambient), h_m_5(i)); % mol/m2*s  
    end
   
    T_drop(i) = T_ambient - sol(i, 1); 

end

% Water usage calculations

molar_H2O_evap_rate = H2O_flux_4.*SA_4 + H2O_flux_5.*SA_5; % mol/s

MM_H2O = 18.015*10^(-3); % kg/mol 
rho_water = 1000; % kg/m3

volume_evap_rate = molar_H2O_evap_rate.*MM_H2O./rho_water; % m3/s 
% Convert from m3 to mL 
volume_evap_rate = volume_evap_rate.*1000*3600; % L/hr 

% Integrate using Trapz to get the total volume lost 
total_water_lost = trapz(t./3600, volume_evap_rate); % [L]


disp(['Total Water Evaporated: ', num2str(total_water_lost), ' [L]']); 


% Plotting
figure(1);  
hold on;
plot(t./3600, sol(:,1) - 273.15, '-', 'linewidth', 2.5); % Modeled AIr temperature
% % % Plot ambient and wetbulb temperatures
plot(t./3600, interp1(t_recorded, T_outside, t), '-.', 'linewidth', 2); % ambient
plot(t_recorded(1:(end))./3600, T_inside(1:(end)),  'linewidth', 2.5); % observed Inner chamber 
plot(t./3600, Twb_vecinterp, '-',  'linewidth', 2); % wet bulb
l1 = legend('$T_{1}$','$T_{inf}$', 'Observed', '$T_{wb}$');
set(l1, 'interpreter', 'latex');
xlabel('\textbf{Time [hours]}', 'interpreter', 'latex');
ylabel('\textbf{Temperature [C]}', 'interpreter', 'latex');
set(gca,'FontSize',16,'LineWidth',2,'FontName','Arial');
set(gcf,'Color','w');
figure(1);


%% Error Calculations
T_model = sol(:,1) - 273.15; % [Celsius]

T_interp = interp1(t, sol(:,1) - 273.15, t_recorded); % modeled data

A = isnan(T_interp);
nan_index_array = []; 
k = 1; 

for i = 1:length(T_interp)
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


T_interp = rmmissing(T_interp); % modeled data 
T_inside = rmmissing(T_inside); % observed 


% Calculation of Error between model and 
T_obs_mean = mean(T_inside); % experimental 
SS_total = sum((T_inside - T_obs_mean).^2); % experimental compared to experimental mean 
SS_res = sum((T_inside -  T_interp).^2);

T_model_mean = mean(sol(:, 1) - 273.15);

NMSE = (SS_res/(T_model_mean*T_obs_mean))/(length(sol(:,1))); % normalized mean squared error

MAE = sum(abs(T_interp - T_inside))/length(T_interp); % Mean absolute error 

RMSE = (SS_res/(length(T_inside)))^(1/2); % root mean squared error

toc







