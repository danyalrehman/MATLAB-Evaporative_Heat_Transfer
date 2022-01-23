function h_vec = ConvectionModel(speeds, T_vec, t, parameters, ambient_vars)

% This function is called for a given time "t" to calculate the convection
% coefficient 


h_vec = zeros(3, 1);
 
% Establish variables correctly depending on whether vegetables are present
V_veges = parameters(2, 3);

% Row 1 is h1 (inner chamber)
% Row 2 is h4
% Row 3 is h5 

% t is the current time step 

% T_vec is the vector of temperatures at time t 
% Unpack temperatures (Kelvin)

T_1 = T_vec(1);

if V_veges == 0
    T_5 = T_vec(5);
    T_4 = T_vec(4);
    T_2 = T_vec(2); 
else
    T_5 = T_vec(6);
    T_4 = T_vec(5);
    T_2 = T_vec(3);  
end

T_ambient = ambient_vars(1);

mu_inf = ambient_vars(3);
k_inf = ambient_vars(4);
Cp_inf = ambient_vars(5);
rho_inf = ambient_vars(6);

D_waterAir_4 = D_wa(T_4);
D_waterAir_5 = D_wa(T_5);


%% Parameters 
% Unpack or calculate all the parameters 
Cp_vec = parameters(:, 2); % Cp
Cp_1 = Cp_sat_air(T_1);

% Densities 
rho_vec = parameters(:,1); % Density 
rho_1 = rho_sat_air(T_1);


% Thermal Conductivity 
k_vec = parameters(:, 6); % W/mK
k1 = k_sat_air(T_1);

% Viscosity 
mu_1b = mu_sat_air(T_1);

% Radii and Heights 

r_vec = parameters(:, 4);
r4 = r_vec(5);

H_vec = parameters(:, 5);
H1 = H_vec(1);
H4 = H_vec(5);


inner_speed = speeds(1);
outer_speed = speeds(2); 

%% Correlations 


g = 9.8; % m/s^2 
alpha_1b = k1/(rho_1*Cp_1); % Thermal Diffusivity of the air in 1b m2/s 
nu_1b = mu_1b/rho_1; % Kinematic Viscosity of Air in section 1b [m^2/s]
nu_inf = mu_inf/rho_inf; % m2/s 
alpha_inf = k_inf/(rho_inf*Cp_inf);
B_outside = 1/T_ambient; % Assume Air behaves as an Ideal Gas 

Pr_1b = nu_1b/alpha_1b; 
Pr_inf = nu_inf/alpha_inf;
n = 7/2; % Nusselt combined forced and free parameter for cylinders 

L5 = H_vec(6); % thickness of lid/cloth 

L_char_5 = 2*r4;

if outer_speed == 0 % Free Convection only 
    Ra_5 = g*B_outside*(T_5 - T_ambient)*(L_char_5^3)/(nu_inf*alpha_inf);
    Ra_5 = abs(Ra_5);
    
    Nu_5 = 0.52*Ra_5^(1/5); 

    % For the Sides (Region 4/outside boundary)  

    Ra_4 = g*B_outside*(T_4 - T_ambient)*((2*H4 + L5)^3)/(nu_inf*alpha_inf);
    Ra_4 = abs(Ra_4);

    if Ra_4 <= 10^9
        Nu_4 = 0.68 + (0.670*Ra_4^(1/4))/(1 + (0.492/Pr_inf)^(9/16))^(4/9); 
    else
        Nu_4 = (0.825 + (0.387*Ra_4^(1/6))/(1 + ((0.492/Pr_inf)^(9/16)))^(8/27))^2;
    end

else % Forced Convection and free convection 
    % Lid (5)
    Re_4 = rho_inf*outer_speed*2*H4/(mu_inf);
    Re_5 = rho_inf*outer_speed*L_char_5/(mu_inf);
    % Schmidt Number 
    Sc_amb_4 = mu_inf/(rho_inf*D_waterAir_4);
    Sc_amb_5 = mu_inf/(rho_inf*D_waterAir_5);
    if Re_5 < 5*10^5 % Laminar 
        Nu_5_Forced = 0.664*(Re_5^(1/2))*Pr_inf^(1/3);
        
        % Mass Transfer 
        
        % Sherwood Number
        Sh_amb_4 = 0.664*(Re_4^(1/2))*Sc_amb_4^(0.3);
        Sh_amb_5 = 0.664*(Re_5^(1/2))*Sc_amb_5^(0.3);
    
    else  % Turbulent
        Sh_amb_4 = 0.037*(Re_4^(0.8))*Sc_amb_4^(0.3);
        Sh_amb_5 = 0.037*(Re_5^(0.8))*Sc_amb_5^(0.3);
        Nu_5_Forced = (Pr_inf^(1/3))*(0.037*(Re_5^(4/5)) - 871); 
    end
    % Wall (4)
    Nu_4_Forced = 0.3 + ((0.62.*(Re_4^(1/2))*(Pr_inf^(1/3))/(1 + ...
        (0.4/Pr_inf)^(2/3))^(1/4)))*(1 + (Re_4/(282000))^(5/8))^(4/5); 
    
    % Mass Transfer coefficient 
    h_m_4 = Sh_amb_4*D_waterAir_4/(2*H4 + L5);
    h_m_5 = Sh_amb_5*D_waterAir_5/(L_char_5);
    
    % Free Convection (Lid)
    Ra_5 = g*B_outside*(T_5 - T_ambient)*(L_char_5^3)/(nu_inf*alpha_inf);
    Ra_5 = abs(Ra_5);
    Nu_5_Free = 0.52*Ra_5^(1/5);

    % For the Sides (Region 4/outside boundary)  

    Ra_4 = g*B_outside*(T_4 - T_ambient)*((2*H4 + L5)^3)/(nu_inf*alpha_inf);
    Ra_4 = abs(Ra_4);

    if Ra_4 <= 10^9
        Nu_4_Free = 0.68 + (0.670*Ra_4^(1/4))/(1 + (0.492/Pr_inf)^(9/16))^(4/9); 
    else
        Nu_4_Free = (0.825 + (0.387*Ra_4^(1/6))/(1 + ((0.492/Pr_inf)^(9/16)))^(8/27))^2;
    end
    
    % Get the combined nusselt number
    
    Nu_4 = (Nu_4_Forced^n + Nu_4_Free^n)^(1/n); 
    Nu_5 = (Nu_5_Forced^n + Nu_5_Free^n)^(1/n);  
    
end

h5 = Nu_5*k_inf/(L_char_5);
h4 = Nu_4*k_inf/(2*H4 + L5);

%% Inner Chamber 

 % For the inner chamber 

if inner_speed == 0 % Free Convection
    B_1 = 1/(T_1);
    Ra_1 = g*B_1*(T_2 - T_1)*((2*H1)^3)/(nu_1b*alpha_1b);
    Ra_1 = abs(Ra_1);
    Nu_1 = 0.55*Ra_1^(1/4);
    h1 = Nu_1*k1/(2*H1);
else   
    % Forced 
    
    Re_1 = rho_1*inner_speed*2*r4/(mu_1b);
    % Forced Convection
    Nu_1_Forced = 0.3 + ((0.62.*(Re_1^(1/2))*(Pr_1b^(1/3))/(1 + ...
        (0.4/Pr_1b)^(2/3))^(1/4)))*(1 + (Re_1/(282000))^(5/8))^(4/5); 
    
    % Free Convection
    B_1 = 1/(T_1);
    Ra_1 = g*B_1*(T_2 - T_1)*((2*H1)^3)/(nu_1b*alpha_1b);
    Ra_1 = abs(Ra_1);
    Nu_1_Free = 0.55*Ra_1^(1/4);
    
    Nu_1 = (Nu_1_Free^n + Nu_1_Forced^n)^(1/n);
    
    h1 = Nu_1*k1/(2*H1);
    
end

h_vec(1) = h1; 
h_vec(2) = h4;
h_vec(3) = h5;

end