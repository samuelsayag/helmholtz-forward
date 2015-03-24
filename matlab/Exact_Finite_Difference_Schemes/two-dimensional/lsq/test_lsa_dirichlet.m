%==========================================================================
% Test of the algorithm that will optimize the theta
%==========================================================================

% close all; clear all; clc;
addpath(genpath('..\..\..\..\matlab'));

% definition of the dirichlet function with the search theta
s_theta = pi/4;
p.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    s_theta, i * params.h, (j-1) * params.h);
p.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    s_theta, (i-1) * params.h, j * params.h);
p.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    s_theta, i * params.h, (j+1) * params.h);
p.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    s_theta, (i+1) * params.h, j * params.h);


% parameter of a helmholtz simulation
p.h = [0.02];
p.k = 300 * sqrt(2);
p.a = 0;
p.b = 1;
p.d = 1;
p.c = 0;
p.m = (p.b - p.a)./p.h;
p.n = (p.d - p.c)./p.h;
p.interior = 'new';
% p.boundary = 'new';

% STEP 2.0 - ????? how to take the initial value of theta that we want to use for
% sommerfeld. in the doubt 
% p.theta = pi/3;

% building of the function pointer
x = linspace(1, p.m, p.m) * p.h;
y = linspace(p.n, 1, p.n) * p.h;
[X,Y] = meshgrid(x,y);
af = @(t) analytic_sol_2D(p.k, t, X, Y);

% parameter of the non linear least square algorithm
lsp.theta1 = 0;
lsp.theta2 = pi/2;
lsp.na     = 50;
lsp.pr     = 1e-3;

res = lsa_3(p, af, lsp);