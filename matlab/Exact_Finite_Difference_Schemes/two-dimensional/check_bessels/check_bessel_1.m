%==========================================================================
% script that shows that the computed function of bessel by means of
% integral is giving exactly the same result as the function of besselj
% given natively by the matlab software.
%==========================================================================

close all; clear all; clc;
addpath(genpath('..\..\..\..\matlab'));
pause on;

% generic parameters of the simulations
k = [10, 20, 30, 50, 100, 150];
d_k = size(k, 2);
sim_param.theta = pi/4;

sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.m = 100;
sim_param.n = 100;
sim_param.h = (sim_param.d - sim_param.c)./(sim_param.m - 1);

% dirichlet boundary
% parameters necessary to compute boundary points
sim_param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j-1) * params.h);
sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j+1) * params.h);
sim_param.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i+1) * params.h, j * params.h);

% declaration of solution structures
sols = {};
params = {};

% simu interior NEW
sim_param.interior = 'new';
sim_param.bessel = @(x) besselj(0, x);
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.bessel = @(x) bessel_integral(x, 0, pi);
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];


% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
c_k = mat2cell(k', ones(1, d_k));
x = linspace(1, sim_param.m, sim_param.m) * sim_param.h;
y = linspace(sim_param.n, 1, sim_param.n) * sim_param.h;
[X,Y] = meshgrid(x,y);

% calculate the error for each simulation
error = cell(size(sols));
analytic = cell(size(sols));
fa = @(t) analytic_sol_2D(t, sim_param.theta, X, Y);
fe = @(a, s) norm((a - s), Inf );
for j = 1:size(error,2)
    analytic(:,j) = cellfun(fa, c_k, 'UniformOutput', false);
    error(:,j) = cellfun( fe, analytic(:,j), sols(:,j), 'UniformOutput', false );
end

% % ONE COLUMN RESULT
% % just the central scheme
title1 = {'' 'error' 'error' };
title2 = {'k' 'NC-NB - bessel (matlab)' 'NC-NB - bessel (integral)'};
res_tab = [title1; title2];
res_tab = [res_tab; c_k  error ];
res_tab