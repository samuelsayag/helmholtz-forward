%--------------------------------------------------------------------------
% Build and plot analytic solution of Helmholtz equation from article:
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION
% Yogi Erlangga, and Eli Turkel 
%--------------------------------------------------------------------------
% close all; clear all; clc;

% define the parameter of the solution
k = 10;
h = 0.01;

a = 0; b = 1;
c = -1/2; d = 1/2;

m = (d-c)/h + 1;
n = (b-a)/h + 1;

x_v = linspace(a, b, m);
y_v = linspace(d, c, n);
[X,Y] = meshgrid(x_v,y_v);

Z = helm_sol1(X, Y, k);

figure(2)
subplot(1, 2, 1);
mesh(X, Y, real(Z));
title 'Real part'
subplot(1, 2, 2);
mesh(X, Y, imag(Z));
title 'Img part'