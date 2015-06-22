%==========================================================================
% Test of the algorithm that will optimize the theta
%==========================================================================

close all; clear variables; clc;
addpath(genpath('..\..\..\..\matlab'));

% definition of the dirichlet function with the search theta
p.theta = pi/4;
p.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, (j-2) * params.h);
p.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-2) * params.h, (j-1) * params.h);
p.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, (j) * params.h);
p.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i) * params.h, (j-1) * params.h);


% parameter of a helmholtz simulation
p.h = [0.02];
p.k = 60 * sqrt(2);
p.a = 0;
p.b = 1;
p.d = 1;
p.c = 0;
p.m = (p.b - p.a)./p.h+1;
p.n = (p.d - p.c)./p.h+1;
p.interior = 'new';
% p.boundary = 'new';

% STEP 2.0 - ????? how to take the initial value of theta that we want to use for
% sommerfeld. in the doubt 
% p.theta = pi/3;

% building of the function pointer
% analytic solution
x = linspace(p.a, p.b, p.m);
y = linspace(p.d, p.c, p.n);
[X,Y] = meshgrid(x,y);
af = @(t) analytic_sol_2D(p.k, t, X, Y);

% parameter of the non linear least square algorithm
lsp.theta1 = pi/16;
lsp.theta2 = 7*pi/16;
lsp.na     = 10;
lsp.pr     = 1e-6;

res = lsa_4(p, af, lsp);