function Emiss_clear = E_clear(T_dp)

a = T_dp/100; 

Emiss_clear = 0.711 + 0.56*(a + 0.73*a^2);



end