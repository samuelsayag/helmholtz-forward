% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

% clear all;
close all; clear variables; clc;

% modeled solution
theor = @(x, y, k) helm_sol1_2D( x, y, k );

% basic parameter of the simulation
param.k = 30;
% definition of the place
param.a = 0; 
param.b = 1;
param.c = -1/2; 
param.d = 1/2;
size = 50;
param.m = size;
param.n = size;
param.h = (param.d - param.c)/(param.m-1);
% dirichlet function
param.dirichlet = @(x,y) theor( x, y, param.k );
scheme = Ord6thHelmholtz2D(param.k, param.h);

% define the solver
solver = @(A, b) A\b;

param
ps = ProblemSolver(param, scheme, solver);
[ A, b, sol ] = ps.solve();

error = ErrorHandler( param, sol );
error

axis_scale = [param.a, param.b, param.c, param.d, -1, 1];
compare_graphs( param,  sol, error, axis_scale, 1);
