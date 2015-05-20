function [ area ] = bessel_integral( x, a, b, nb_pt )
% BESSEL_INTEGRAL 
% Compute a Bessel function from the integral proposed in the article
% EXACT FINITE DIFFERENCE SCHEMES FOR SOLVING HELMHOLTZ EQUATION AT ANY 
% WAVENUMBER, YAU SHU WONG AND GUANGRUI LI, 2011.
% 
% param:
%   x: the abscissa
%   a: (optional) the lower bound of the integral
%   b: (optional) the higher bound of the integral
% return:
%   area: the integral

pt = 1e3;

if nargin < 1
    error('not enough parameters (at least one)')
elseif nargin < 2
    a = 0; b = pi; nb_pt=pt;
elseif nargin < 3
    b = pi; nb_pt=pt;
elseif nargin < 4
    nb_pt=pt;
end

if a > b 
    error('a may not be superior to b in the computation of the integral')
end

fun = @(theta) bessel_exact_theta(x, theta);
X = linspace(a, b, nb_pt);
Y = fun(X);
area = (1/(b-a)) * trapz(X,Y);

end

