function Emiss_clear = E_clear(T_dp)

% Berdahl and Martin (1984) empirical correlation 
% https://scholarworks.uark.edu/cgi/viewcontent.cgi?referer=https://www.google.com/&httpsredir=1&article=2230&context=etd



a = T_dp/100; 

Emiss_clear = 0.711 + 0.56*(a + 0.73*a^2);



end