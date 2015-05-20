%--------------------------------------------------------------------------
% Build and plot analytic solution of Helmholtz equation from article:
% EXACT FINITE DIFFERENCE SCHEMES FOR SOLVING HELMHOLTZ EQUATION AT ANY 
% WAVENUMBER, YAU SHU WONG AND GUANGRUI LI, 2011.
%--------------------------------------------------------------------------
%close all; clear all; clc;

% define the parameter of the solution
k = 30 * sqrt(2);
theta = pi/4;
% define the region we want to have the solution on
h = 0.01;
a = 0;
b = 1;
c = 1;
m = (b-a)/h + 1;
n = (c-a)/h + 1;
x = linspace(0,1, m);
y = linspace(1,0, n);
[X,Y] = meshgrid(x,y);

Z = analytic_sol_2D(k, theta, X, Y);

figure(1)
subplot(1, 2, 1);
mesh(X, Y, real(Z));
title 'Real part'
subplot(1, 2, 2);
mesh(X, Y, imag(Z));
title 'Img part'