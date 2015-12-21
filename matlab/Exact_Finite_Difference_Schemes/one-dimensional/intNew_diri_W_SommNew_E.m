function [A, b, x, param] = intNew_diri_W_SommNew_E()
% std_diriW_diriE
% Summary of this function goes here
%   Detailed explanation goes here
%close all; clear all; clc;

% x = [];
% 
% 
% parameters necessary to compute interior points
% params.k = 10;
% params.h = 1e-2;
% params.h = 0.1;
% a = 0;
% b = 1;
% 
% forcingDirichlet = true;
% 
% params.m = (b-a)./params.h;
% 
% params.interior = 'new';
% params.boundary = 'sommerfeld_new';
% params.dirichlet.W = @(params, A, b, i) 1;
% 
% 
% params
% 
% [ func_scheme, params ] = helmholtz_1D_scheme_factory( params );
% 
% create the matrix of finite difference
% [A, b] = build_one_dimensional_problem2(params, func_scheme);
% 
% % debug
% full(A)
% b
% 
% ------------------- solve the system -----------------------------
% tstart = tic;
% 
% tol = 1e-6;
% [x,flag,relres] = gmres(A, b, [], tol);        
% [x,flag,relres] = gmres(A, b);
% x = inv(A) * b;
% x = [1;x];
% 
% telapsed = toc(tstart);
% ------------------- display some result -----------------------------
% 
% telapsed

close all; clear all; clc;
%------------------------- parameter of the simulation --------------------
% define the domain of the BVP
a0 = 0;
a1 = 1;
% parameters of the BVP
param.k = 10;
param.h = 1e-2;
param.m = (a1-a0)./param.h + 1;
param.interior = 'new';
param.boundary = 'sommerfeld_new';
% Dirichlet boundary on the West side, the factory deduce it has to use
% Sommerfeld on the East side. 
% Here the parameter are not used because we use a constant value but the
% user may potentially use this parameter for a Dirichlet condition.
param.dirichlet.W = @(params, A, b, i) analytic_sol_1D(param.k, a0-param.h);
%param.dirichlet.E = @(params, A, b, i) analytic_sol_1D(param.k, a1+param.h);

%------------------------- use of the factory -----------------------------
% Instantiate the scheme function from the given parameters
[ func_scheme, param ] = helmholtz_1D_scheme_factory( param );

%------------------------- building the matrix ----------------------------
% create the matrix of finite difference
[A, b] = build_one_dimensional_problem2(param, func_scheme);

%------------------------- solving the system -----------------------------
x = A\b;

%------------------------- checking output --------------------------------
% building the close problem value on the domain
d = linspace(a0, a1, param.m);
sol = analytic_sol_1D(param.k, d);

% show a fraction of the vector
tit = {'close sol' 'computed sol' '(last 10 values)'};
disp(tit);
disp([sol(92:end) x(92:end)]);

% error on the whole vector
normInf = norm(sol - x, inf);
disp('Infinite norm');
disp(normInf);