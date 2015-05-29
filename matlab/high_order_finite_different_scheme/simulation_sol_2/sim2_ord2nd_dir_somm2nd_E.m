% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

clear variables; close all; clc;

% modeled solution
theor = @(x, y, k, theta) helm_sol2_2D( k, theta, x, y);
theta = pi/4;

% basic parameter of the simulation
param.k = 10;
param.h = 0.02;
% definition of the area we simulate in it
param.a = 0; 
param.b = 1;
param.c = 0; 
param.d = 1;
param.m = (param.d - param.c)/param.h + 1;
param.n = (param.b - param.a)/param.h + 1;

% boundary condition
param.dirichlet = @(x,y) theor( x, y, param.k , theta);
scheme = Ord2ndHelmholtz2D(param.k, param.h);
param.east = 'sommerfeld';
beta = - param.k * cos(theta);
sommerfeld = Ord2ndSommerfeld2D( param.h, beta);

% define the solver
solver = @(A, b) A\b;

param
ps = ProblemSolver(param, scheme, solver, sommerfeld);
[ A, b, sol ] = ps.solve();

error = ErrorHandler( param, sol );
error

axis_scale = [param.a, param.b, param.c, param.d, -1, 1];
compare_graphs( param,  sol, error, axis_scale, 1);