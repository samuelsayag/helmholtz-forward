% SIMULATIO OF THE HELMHOLTZ EQUATION
% operate a simulation of the 5*5 region
% Dirichlet bound on the south and west side
% Sommerfeld on the north and east side

%close all; clear all; clc;

params.k = 10;
params.h = 0.25;
% parameters of the region
a = 0; 
b = 1;
c = 0;
d = 1;
params.m = (b-a)/params.h;
params.n = (d-c)/params.h;
% parameters necessary to compute boundary points
% params.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, i * params.h, (j-1) * params.h);
params.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
params.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j+1) * params.h);
% params.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i+1) * params.h, j * params.h);

params.interior = 'new';
params.theta = pi/4;
params.boundary = 'new';
% params.bessel = @(x) bessel_exact_theta(x, params.theta);
% params.bessel = @(x) besselj(0, x);

% create the scheme function
[ func_scheme, params ] = helmholtz_2D_scheme_factory( params );
% create the matrix of finite difference
[A, b] = build_two_dimensional_problem2(params, func_scheme);

% full(A)
% full(b)

% x = gmres(A,b');
x = A\b;
sol = transpose(reshape(x, params.m, params.n));

% analytic solution
x = linspace(1, params.m, params.m) * params.h;
y = linspace(params.n, 1, params.n) * params.h;
[X,Y] = meshgrid(x,y);

% build experimental sol with its boundary
ana = analytic_sol_2D(params.k, params.theta, X, Y);

display('norm inf ');
norm((ana-sol), inf)

figure(1)
subplot(1, 2, 1);
mesh(X, Y, real(ana));
title('analytic solution')
subplot(1, 2, 2);
mesh(X, Y, real(sol));
title('computed solution')


