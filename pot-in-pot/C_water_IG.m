function C_H2O = C_water_IG(T, RH)

% T in Kelvin; RH as a decimal 

P_H2O = RH.*Saturation_pressure(T);

R = 8.314; % J/molK 

% Ideal Gas Law 
C_H2O = P_H2O./(R.*T);


end