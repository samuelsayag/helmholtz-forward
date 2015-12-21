% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

clear variables; close all; clc;

% modeled solution
theta = pi/4;
theor = @(x, y, k, theta) helm_sol2_2D( k, theta, x, y);

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

% Dirichlet function and Sommerfeld boundary
param.dirichlet = @(x,y) theor( x, y, param.k , theta);
scheme = ExactScheme2D(param.k, param.h, theta);
param.north = 'sommerfeld';
param.west = 'sommerfeld';
beta = - param.k;
sommerfeld_scheme = ExactSommerfeld2D( param.h, beta, theta, scheme);

% define the parallelisation mode
param.parall = 'parall_2';

% define the solver
solver = @(A, b) A\b;

display(param)
ps = ProblemSolver(param, scheme, solver, sommerfeld_scheme);
[ A, b, sol ] = ps.solve();

error = ErrorHandler( param, sol );
display(error);

axis_scale = [param.a, param.b, param.c, param.d, -1, 1];
compare_graphs( param,  sol, error, axis_scale, 1 );