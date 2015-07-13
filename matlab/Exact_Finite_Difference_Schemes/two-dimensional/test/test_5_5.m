close all; clear variables; clc;
addpath(genpath('..\..\..\matlab'));
pause on;

% generic parameters of the simulations
% sim_param.h = 1e-2;
sim_param.h = 0.2;
sim_param.k = [5];
% sim_param.k = [30 * sqrt(2)];

sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.m = (sim_param.b - sim_param.a)./sim_param.h +1;
sim_param.n = (sim_param.d - sim_param.c)./sim_param.h +1;
sim_param.theta = pi/4;
% dirichlet boundary
% parameters necessary to compute boundary points
sim_param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, (j-2) * params.h);
sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-2) * params.h, (j-1) * params.h);
% sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i-1) * params.h, (j) * params.h);
% sim_param.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i) * params.h, (j-1) * params.h);

sim_param.interior = 'new';
sim_param.boundary = 'new';

sim_param

% build the simulation function
[ func_scheme, sim_param ] = helmholtz_2D_scheme_factory( sim_param );
% generate the matrix and vector of the problem
[At, bt] = build_two_dimensional_problem2(sim_param, func_scheme);
% compute the solution of the equation
solt = At\bt;
solt = transpose(reshape(solt, sim_param.m, sim_param.n));

x = linspace(sim_param.a, sim_param.b, sim_param.m);
y = linspace(sim_param.d, sim_param.c, sim_param.n);
[X,Y] = meshgrid(x,y);

analytic = analytic_sol_2D(sim_param.k, sim_param.theta, X, Y);
error = norm((analytic - solt), Inf );

error