% simulation with variable k ord6th scheme.

clear variables; close all; clc;

% parameter for the component (k, beta, sol) of the closed solution
a = 5;
b = 4;
c = -10;
[ pSol, k, beta ] = analytic_sol_3( a, b, c );

% parameters for the grid on it we model the solution
% we model the problem on the square: 
% x belong [a,b],  y belong [c,d]
param.a = 0;  
param.b = pi;
param.c = 0; 
param.d = pi;
param.n = 300; % number of column of the grid
param.m = 300; % number of line of the grid
% we deduce the step of the grid from this data
param.h = (param.b - param.a)/ (param.n - 1);
% dirichlet function
param.dirichlet = pSol;


% -------------------------------------------------------------------------
% computation of the derivative from the matrix of k(x,y)
% we want to demonstrate the use of the matrix derivative function
% -------------------------------------------------------------------------
dx = linspace(param.a - param.h, param.b + param.h, param.n+2);
dy = linspace(param.d + param.h, param.c - param.h, param.m+2);
[dX, dY] = meshgrid(dx,dy);
k_mat = k(dX, dY);
% we may suppose now we dispose of an empirical matrix k_mat of dimension
% [m+1 x n+1] and we will build the serie of its derivative.
[k2, K2x, K2y, K2xx, K2yy] = ...
    Ord6thVarKHelmholtz2D.build_matricial_derivative(k_mat, param);
scheme = Ord6thVarKHelmholtz2D( param.h, k2, K2x, K2y, K2xx, K2yy);

% define the solver
solver = @(A, b) A\b;

param
ps = ProblemSolver(param, scheme, solver);
[ A, b, sol ] = ps.solve();

error = ErrorHandler( param, sol );
error

axis_scale = [param.a, param.b, param.c, param.d, -5, 5];
compare_graphs( param,  sol, error, axis_scale, 1);