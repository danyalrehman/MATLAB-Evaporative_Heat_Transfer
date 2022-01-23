 function [Pws] = Saturation_pressure(T) %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6
        % Input in Kelvin
        Pws=exp(-(5.8002206e3)./T+1.3914993+-(4.8640239e-2).*T+(4.1764768e-5).*(T.^2)-(1.4452093e-8).*(T.^3)+6.5459673.*log(T)); %in Pa valid for 0 to 200C
        % in Pa
 end