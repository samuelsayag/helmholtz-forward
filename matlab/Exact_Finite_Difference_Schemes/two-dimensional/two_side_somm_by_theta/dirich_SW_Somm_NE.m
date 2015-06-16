%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\..\..\matlab'));

% generic parameters of the simulations
sim_param.h = [0.02];
% k = [25];
% k = [150, 100, 70];
% k = sqrt(2)* [30, 25, 20];
k = sqrt(2)* [30, 25, 20, 15, 10, 5];
d_k = size(k,2);


sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.m = (sim_param.b - sim_param.a)./sim_param.h + 1;
sim_param.n = (sim_param.d - sim_param.c)./sim_param.h + 1;

sim_param.theta = pi/4;

% dirichlet boundary
% parameters necessary to compute boundary points
sim_param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, (j-2) * params.h);
sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-2) * params.h, (j-1) * params.h);
% params.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i-1) * params.h, (j) * params.h);
% params.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i) * params.h, (j-1) * params.h);

% declaration of solution structures
sols = {};
params = {};

% simu interior std
sim_param.interior = 'std';
sim_param.boundary = 'new';
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'new';
sim_param.bessel = @(x) bessel_integral(x, 0, pi);
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'new';
sim_param.bessel = @(x) bessel_integral(x, pi/8, 3 * pi/8);
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'new';
sim_param.bessel = @(x) bessel_exact_theta(x, sim_param.theta);
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
c_k = mat2cell(k', ones(1, d_k));
x = linspace(sim_param.a, sim_param.b, sim_param.m);
y = linspace(sim_param.d, sim_param.c, sim_param.n);
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

% RESULT FOR STANDARD AND NEW ('std', 'new')
% just the central scheme
title1 = {'' 'error' 'error' 'error' 'error'};
title2 = { 'k' 'SFD' 'NFD-J0[0,pi]' 'NFD-J0[pi/8,3pi/8]'  'NFD - exact theta'};
res_tab = [title1;title2];
res_tab = [res_tab; c_k error ];
res_tab

