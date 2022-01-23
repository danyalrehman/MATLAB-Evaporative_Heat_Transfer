function D_AB = D_wa(T)

% This function calculates the diffusivity of water in air as a function of
% temperature, inputted in Kelvin

D_AB = -2.755E-6 + 4.479E-8*T + 1.656E-10*T^2;  % m^2/s

end