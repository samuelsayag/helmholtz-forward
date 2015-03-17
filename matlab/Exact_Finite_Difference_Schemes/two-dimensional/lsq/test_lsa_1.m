%==========================================================================
% Test of the algorithm that will optimize the theta
%==========================================================================

% close all; clear all; clc;
addpath(genpath('..\..\..\..\matlab'));

% parameter of a helmholtz simulation
p.h = [0.2];
p.k = 15 * sqrt(2);
p.a = 0;
p.b = 1;
p.d = 1;
p.c = 0;
p.m = (p.b - p.a)./p.h;
p.n = (p.d - p.c)./p.h;
p.interior = 'new';
p.boundary = 'new';
% dirichlet boundary
% parameters necessary to compute boundary points
p.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j-1) * params.h);
p.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
% sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, i * params.h, (j+1) * params.h);
% sim_param.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i+1) * params.h, j * params.h);

% building of the function pointer
x = linspace(1, p.m, p.m) * p.h;
y = linspace(p.n, 1, p.n) * p.h;
[X,Y] = meshgrid(x,y);
af = @(t) analytic_sol_2D(p.k, t, X, Y);

% parameter of the non linear least square algorithm
lsp.theta1 = 0;
lsp.theta2 = pi/2;
lsp.na     = 6;
lsp.pr     = 1e-3;

res = lsa_1(p, af, lsp);