%--------------------------------------------------------------------------
% Build a Bessel function from the integral proposed in the article
% EXACT FINITE DIFFERENCE SCHEMES FOR SOLVING HELMHOLTZ EQUATION AT ANY 
% WAVENUMBER, YAU SHU WONG AND GUANGRUI LI, 2011.
%--------------------------------------------------------------------------
close all; clear all; clc;

h = 0.02;
k = sqrt(2) * [30 25 20 15 10 5];
theta = pi/4;

r = cell(8, 5);
r(1,:) = {'', '', 'J0(kh)', 'J0(kh)', 'J0(kh)'};
r(2,:) = {'kh', 'k', 'sum [0,pi]', 'exact theta', 'matlab'};

coef = 1;

for j=1:size(k,2)
    i = j + 2;
    kh = k(j) * h;
    r{i, 1} = kh;
    r{i, 2} = k(j);
    r{i, 3} = coef * bessel_integral(kh);
    r{i, 4} = coef * bessel_exact_theta(kh, theta);
    r{i, 5} = coef * besselj(0, kh);       
end

r
