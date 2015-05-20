function [ res ] = bessel_exact_theta( x, theta )
%BESSEL_EXACT_THETA 
% Compute the exact theta (at one point) Bessel function describe in
% EXACT FINITE DIFFERENCE SCHEMES FOR SOLVING HELMHOLTZ EQUATION AT ANY 
% WAVENUMBER, YAU SHU WONG AND GUANGRUI LI, 2011.

res = cos(x * sin(theta));

end

