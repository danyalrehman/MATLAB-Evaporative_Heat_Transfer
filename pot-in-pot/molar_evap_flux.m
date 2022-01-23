function N_evap = molar_evap_flux(C_H2O_n, C_H2O_inf, h_m)


% T_inf and T_s in Kelvin
 


% Concentration of water at surface (Assumed saturated thin film) 
% C_H2O_s = C_water_IG(T_s, 1);


% D_AB at 310K ~ 2.76e-05
% Le = 0.8699; % Lewis number at 310K 

% Combined Reynold's Analogy 
% h_m = h/(rho_a*Cp_a*Le^(2/3)); % [m/s]

% delta = 0.0083; % [m] 
% k = 6E-12; % permeability % 
% P_n = Saturation_pressure(T_s); % Pa 
% P_inf = RH_inf*Saturation_pressure(T_inf); % Pa 
% MM_H2O = 18.015*10^(-3); % kg/mol 
% 
% 
% N_evap = (k/delta)*(P_n - P_inf)/MM_H2O; % [mol/m2*s]

N_evap = h_m.*(C_H2O_n - C_H2O_inf); % mol/m2*s 


end