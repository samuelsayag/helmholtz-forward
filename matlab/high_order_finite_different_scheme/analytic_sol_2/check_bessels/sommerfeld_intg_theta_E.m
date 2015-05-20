%==========================================================================
% simulation with a fixed theta and a redefine integral bessel function
% redefine for each simulation (each different theta)
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\..\..\..\matlab'));
pause on;

% generic parameters of the simulations
sim_param.h = [0.02];
angle_div = 23;
theta = pi/2 * 1/angle_div * linspace(0, angle_div, angle_div + 1);
% theta =   0.2732 + abs(0.4098 - 0.2732) * 1/angle_div * linspace(0, angle_div, angle_div + 1);
d_theta =size(theta,2);

sim_param.k = 30 * sqrt(2);
sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.m = (sim_param.b - sim_param.a)./sim_param.h;
sim_param.n = (sim_param.d - sim_param.c)./sim_param.h;
% dirichlet boundary
% parameters necessary to compute boundary points
sim_param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j-1) * params.h);
sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
% sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, i * params.h, (j+1) * params.h);

% declaration of solution structures
sols = {};
params = {};

tic

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'new';

% sim_param.bessel = @(x) bessel_integral(x, 0, pi/2);
% [ sol, param ] = simulation_theta_2D( theta, sim_param );

[ sol, param ] = simulation_theta_2D( theta, sim_param );

% [ sol, param ] = simulation_intg_theta_2D( theta, sim_param );

sols = [sols, sol];
params = [params, param];

% time elapsed
elapsed = toc;
elapsed

% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
res_theta = mat2cell(theta', ones(1, d_theta));
x = linspace(1, sim_param.m, sim_param.m) * sim_param.h;
y = linspace(sim_param.n, 1, sim_param.n) * sim_param.h;
[X,Y] = meshgrid(x,y);

% calculate the error for each simulation
error = cell(size(sols));
analytic = cell(size(sols));
fa = @(t) analytic_sol_2D(sim_param.k, t, X, Y);
% fa = @(t) analytic_sol_2D(sim_param.k, pi/4, X, Y);
fe = @(a, s) norm((a - s), Inf );
for j = 1:size(error,2)
    analytic(:,j) = cellfun(fa, res_theta, 'UniformOutput', false);
    error(:,j) = cellfun( fe, analytic(:,j), sols(:,j), 'UniformOutput', false );
end

% % ONE COLUMN RESULT
% % just the central scheme
title1 = {'' '' 'error' };
title2 = {'theta' 'cos(theta)' 'NFD - NBC'};
res_tab = [title1; title2];
res_cos_theta = mat2cell(cos(theta'), ones(1, d_theta));
res_tab = [res_tab; res_theta res_cos_theta error ];
res_tab

