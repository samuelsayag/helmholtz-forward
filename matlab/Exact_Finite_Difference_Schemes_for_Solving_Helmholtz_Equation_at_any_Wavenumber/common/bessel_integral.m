function [ area ] = bessel_integral( x )
% BESSEL_INTEGRAL 
% Compute a Bessel function from the integral proposed in the article
% EXACT FINITE DIFFERENCE SCHEMES FOR SOLVING HELMHOLTZ EQUATION AT ANY 
% WAVENUMBER, YAU SHU WONG AND GUANGRUI LI, 2011.
% 
% param:
%   x: the abscissa
% return:
%   area: the integral

fun = @(theta) bessel_exact_theta(x, theta);
area = (1/pi) * integral(fun, 0, pi);
% area = integral(fun, 0, pi);

end

