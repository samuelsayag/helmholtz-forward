% SIMULATIO OF THE HELMHOLTZ EQUATION
% operate a simulation of the 5*5 region
% Dirichlet bound on the south and west side
% Sommerfeld on the north and east side

close all; clear all; clc;

params.k = 30 * sqrt(2);
params.h = 0.02;
% parameters of the region
a = 0; 
b = 1;
c = 0;
d = 1;
params.m = (b-a)/params.h;
params.n = (d-c)/params.h;
% parameters necessary to compute boundary points
params.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j-1) * params.h);
params.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
params.interior = 'new';
params.theta = pi/4;
params.boundary = 'new';

% create the scheme function
func_scheme = helmholtz_2D_scheme_factory( params );
% create the matrix of finite difference
[A, b] = build_two_dimensional_problem2(params, func_scheme);

% full(A)
% full(b)

% x = gmres(A,b');
x = A\b';
x_sol_partial = reshape(x, params.m, params.n)';

% analytic solution
x = linspace(0,1, params.m+1);
y = linspace(1,0, params.n+1);
[X,Y] = meshgrid(x,y);
Z = analytic_sol_2D(params.k, params.theta, X, Y);

% build experimental sol with its boundary
x_sol = analytic_sol_2D(params.k, params.theta, X, Y);
x_sol(1:end-1, 2:end) = x_sol_partial(:,:);

figure(1)
subplot(1, 2, 1);
mesh(X, Y, real(Z));
title('analytic solution')
subplot(1, 2, 2);
mesh(X, Y, real(x_sol));
title('computed solution')

display('norm inf = max(sum(abs( Z - ))) : ');
norm(Z-x_sol, inf)
