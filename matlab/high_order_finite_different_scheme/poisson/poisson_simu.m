% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

clear variables; close all; clc;

% definition of the area we simulate in it
param.h = 0.125;%3 pt
% param.h = 0.025;%10 pt

param.a = 0.125; 
param.b = 0.375;
param.c = 0.125; 
param.d = 0.375;
param.m = (param.d - param.c)/param.h + 1;
param.n = (param.b - param.a)/param.h + 1;

% dirichlet function
param.dirichlet = @(x,y) poisson_dirichlet( x, y);
scheme = Poisson2D();

% define the solver
solver = @(A, b) A\b;

param
ps = ProblemSolver(param, scheme, solver);
[ A, b, sol ] = ps.solve();

% graphical representation
x = linspace(param.a,param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid( x, y );
axis_scale = [param.a param.b param.c param.d, 0, 110];

figure(1)
mesh(X, Y, real(sol));
axis(axis_scale)
t1 = sprintf('Heat repartition');
title(t1)