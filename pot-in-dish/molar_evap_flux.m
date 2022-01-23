function N_evap = molar_evap_flux(T_s, T_inf, RH_inf, h, rho_a, Cp_a, Le)


% T_inf and T_s in Kelvin
 

% Concentration of water in ambient air 
C_H2O_inf = C_water_IG(T_inf, RH_inf);

% Concentration of water at surface (Assumed saturated thin film) 
C_H2O_n = C_water_IG(T_s, 1);

% Combined Reynold's Analogy 
h_m = h/(rho_a*Cp_a*Le^(2/3)); % [m/s]



N_evap = h_m.*(C_H2O_n - C_H2O_inf); % mol/m2*s 


end