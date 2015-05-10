%--------------------------------------------------------------------------
% Build and plot analytic solution of Helmholtz equation from article:
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION
% Yogi Erlangga, and Eli Turkel 
%--------------------------------------------------------------------------
close all; clear all; clc;

% define the parameter of the solution
k = 10;
h = 0.1;

a = 0; b = 1;
c = -1/2; d = 1/2;

m = (d-c)/h + 1;
n = (b-a)/h + 1;

x = linspace(0-h,1+h, m+2);
y = linspace(1/2+h, -1/2-h, n+2);
[X,Y] = meshgrid(x,y);

Z = helm_sol1(X, Y, k);

figure(1)
subplot(1, 2, 1);
mesh(X, Y, real(Z));
title 'Real part'
subplot(1, 2, 2);
mesh(X, Y, imag(Z));
title 'Img part'