function h_vec = ConvectionModel(speeds, T_vec, t, parameters, ambient_vars)

% This function is called for a given time "t" to calculate the convection
% coefficient 


h_vec = zeros(5, 1);
 
% Establish variables correctly depending on whether vegetables are present
V_veges = parameters(2, 3);

% t is the current time step 

% T_vec is the vector of temperatures at time t 
% Unpack temperatures (Kelvin)

T_1 = T_vec(1);

if V_veges == 0
    T_5 = T_vec(5);
    T_4 = T_vec(4);
    T_3 = T_vec(3);
    T_2 = T_vec(2); 
else
    T_5 = T_vec(6);
    T_4 = T_vec(5);
    T_3 = T_vec(4);
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
mu_1 = mu_sat_air(T_1);

% Radii and Heights 

r_vec = parameters(:, 4);
r2 = r_vec(3);
r3 = r_vec(4);
r4 = r_vec(5);
t_s = r_vec(6);

H_vec = parameters(:, 5);
H_air = H_vec(1); % Air height
H2 = H_vec(3);
H3 = H_vec(4);
H4 = H_vec(5);


inner_speed = speeds(1);
outer_speed = speeds(2); 

%% Correlations 


g = 9.8; % m/s^2 
alpha_1 = k1/(rho_1*Cp_1); % Thermal Diffusivity of the air in 1b m2/s 
nu_1 = mu_1/rho_1; % Kinematic Viscosity of Air in section 1b [m^2/s]
nu_inf = mu_inf/rho_inf; % m2/s 
alpha_inf = k_inf/(rho_inf*Cp_inf);
B_outside = 1/T_ambient; % Assume Air behaves as an Ideal Gas 

Pr_1 = nu_1/alpha_1; 
Pr_inf = nu_inf/alpha_inf;

n = 7/2; % Nusselt combined forced and free parameter for cylinders 


L5 = H_vec(6); % thickness of lid/cloth 

L_char_5 = 2*r2; % Characteristic length of r2 

if outer_speed == 0 % Free Convection only 
    Ra_5 = g*B_outside*(T_5 - T_ambient)*(L_char_5^3)/(nu_inf*alpha_inf);
    Ra_5 = abs(Ra_5);
    Nu_5 = 0.52*Ra_5^(1/5); 

    % For the Sides (Region 4/outside boundary)  

    Ra_4 = g*B_outside*(T_4 - T_ambient)*(H4^3)/(nu_inf*alpha_inf);
    Ra_4 = abs(Ra_4);

    if Ra_4 <= 10^9
        Nu_4 = 0.68 + (0.670*Ra_4^(1/4))/(1 + (0.492/Pr_inf)^(9/16))^(4/9); 
    else
        Nu_4 = (0.825 + (0.387*Ra_4^(1/6))/(1 + ((0.492/Pr_inf)^(9/16)))^(8/27))^2;
    end

else % Forced Convection and free convection 
   
    Re_2 = rho_inf*outer_speed*(H2 - (H3 - t_s))/(mu_inf);
    Re_3 = rho_inf*outer_speed*(r4 - r3)/(mu_inf);
    Re_4 = rho_inf*outer_speed*H4/(mu_inf);
    Re_5 = rho_inf*outer_speed*L_char_5/(mu_inf);
    % Schmidt Number 
    Sc_amb_4 = mu_inf/(rho_inf*D_waterAir_4);
    Sc_amb_5 = mu_inf/(rho_inf*D_waterAir_5);
    if Re_5 < 5*10^5  
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
    
     % Dish (4) 
    if Re_4 < 2300
       Nu_4_Forced = 0.664*(Pr_inf^(1/3))*Re_4^(1/2);
    else
       Nu_4_Forced = 0.3 + ((0.62.*(Re_4^(1/2))*(Pr_inf^(1/3))/(1 + ...
         (0.4/Pr_inf)^(2/3))^(1/4)))*(1 + (Re_4/(282000))^(5/8))^(4/5); 
    end
    
    
    if Re_3 < 2300
       Nu_3_Forced = 0.664*(Pr_inf^(1/3))*Re_3^(1/2);
    else
       Nu_3_Forced = 0.3 + ((0.62.*(Re_3^(1/2))*(Pr_inf^(1/3))/(1 + ...
         (0.4/Pr_inf)^(2/3))^(1/4)))*(1 + (Re_3/(282000))^(5/8))^(4/5); 
    end
    
    if Re_2 < 2300
       Nu_2_Forced = 0.664*(Pr_inf^(1/3))*Re_2^(1/2); 
    else
       Nu_2_Forced = 0.3 + ((0.62.*(Re_2^(1/2))*(Pr_inf^(1/3))/(1 + ...
         (0.4/Pr_inf)^(2/3))^(1/4)))*(1 + (Re_2/(282000))^(5/8))^(4/5); 
    end
 
    
    % Mass Transfer coefficient 
    h_m_4 = Sh_amb_4*D_waterAir_4/(2*H4 + L5);
    h_m_5 = Sh_amb_5*D_waterAir_5/(L_char_5);
    
    % Free Convection (Lid)
    Ra_5 = g*B_outside*(T_5 - T_ambient)*(L_char_5^3)/(nu_inf*alpha_inf);
    Ra_5 = abs(Ra_5);
    Nu_5_Free = 0.52*Ra_5^(1/5);

    % For the Sides (Region 4/outside boundary)  

    Ra_4 = g*B_outside*(T_4 - T_ambient)*(H4^3)/(nu_inf*alpha_inf);
    Ra_4 = abs(Ra_4);
    
    Ra_3 = g*B_outside*(T_3 - T_ambient)*((r3 - r2)^3)/(nu_inf*alpha_inf);
    Ra_3 = abs(Ra_3);
   
    Ra_2 = g*B_outside*(T_2 - T_ambient)*((H2 - (H3 - t_s))^3)/(nu_inf*alpha_inf);
    Ra_2 = abs(Ra_2);

    if Ra_4 <= 10^9
        Nu_4_Free = 0.68 + (0.670*Ra_4^(1/4))/(1 + (0.492/Pr_inf)^(9/16))^(4/9); 
        Nu_3_Free = 0.68 + (0.670*Ra_3^(1/4))/(1 + (0.492/Pr_inf)^(9/16))^(4/9); 
        Nu_2_Free = 0.68 + (0.670*Ra_2^(1/4))/(1 + (0.492/Pr_inf)^(9/16))^(4/9); 
    else
        Nu_4_Free = (0.825 + (0.387*Ra_4^(1/6))/(1 + ((0.492/Pr_inf)^(9/16)))^(8/27))^2;
        Nu_3_Free = (0.825 + (0.387*Ra_3^(1/6))/(1 + ((0.492/Pr_inf)^(9/16)))^(8/27))^2;
        Nu_2_Free = (0.825 + (0.387*Ra_2^(1/6))/(1 + ((0.492/Pr_inf)^(9/16)))^(8/27))^2;
    end
    
    % Calculate the combined nusselt number
    Nu_2 = (Nu_2_Forced^n + Nu_2_Free^n)^(1/n); 
    Nu_3 = (Nu_3_Forced^n + Nu_3_Free^n)^(1/n); 
    Nu_4 = (Nu_4_Forced^n + Nu_4_Free^n)^(1/n); 
    Nu_5 = (Nu_5_Forced^n + Nu_5_Free^n)^(1/n);  
    
end

h2 = Nu_2*k_inf/(H2 - (H3 - t_s));
h3 = Nu_3*k_inf/(r3 - r2);
h5 = Nu_5*k_inf/(L_char_5);
h4 = Nu_4*k_inf/(H4);

%% Inner Chamber 

 % For the inner chamber 

if inner_speed == 0 % Free Convection
    B_1 = 1/(T_1);
    Ra_1 = g*B_1*(T_2 - T_1)*((H_air)^3)/(nu_1*alpha_1);
    Ra_1 = abs(Ra_1);
    Nu_1 = 0.55*Ra_1^(1/4);
    h1 = Nu_1*k1/(H_air);
else % (Free and Forced)  

    Re_1 = rho_1*inner_speed*H_air/(mu_1);
 
    % Forced Convection
    Nu_1_Forced = 0.3 + ((0.62.*(Re_1^(1/2))*(Pr_1^(1/3))/(1 + ...
        (0.4/Pr_1)^(2/3))^(1/4)))*(1 + (Re_1/(282000))^(5/8))^(4/5); 
    
    % Free Convection
    B_1 = 1/(T_1);
    Ra_1 = g*B_1*(T_2 - T_1)*((H_air)^3)/(nu_1*alpha_1);
    Ra_1 = abs(Ra_1);
    Nu_1_Free = 0.55*Ra_1^(1/4);
    
    Nu_1 = (Nu_1_Free^n + Nu_1_Forced^n)^(1/n);
    
    h1 = Nu_1*k1/(H_air);
    
end

h_vec(1) = h1; 
h_vec(2) = h2;
h_vec(3) = h3;
h_vec(4) = h4;
h_vec(5) = h5;

end