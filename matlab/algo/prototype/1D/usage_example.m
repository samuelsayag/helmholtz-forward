function [ A, b, x, param, close_sol ] = usage( )
%USAGE Summary of this function goes here
%   provide a functional example

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
param.dirichlet.W = @(params,A,b,i) analytic_sol_1D(param.k, a0-param.h);

%------------------------- use of the factory -----------------------------
% Instantiate the scheme function from the given parameters
[ func_scheme, param ] = helmholtz_1D_scheme_factory( param );

%------------------------- building the matrix ----------------------------
% create the matrix of finite difference
tic
[A, b] = build_one_dimensional_problem2(param, func_scheme);
matrix_build_time = toc;

%------------------------- solving the system -----------------------------
x = A\b;
solving_system_time = toc;

%------------------------- checking output --------------------------------
% building the close problem value on the domain
d = linspace(a0, a1, param.m);
close_sol = analytic_sol_1D(param.k, d);

% show a fraction of the vectors of the solution and the computer one
tit = {'close sol' 'computed sol' '(last 10 values)'};
disp(tit);
disp([close_sol(92:end) x(92:end)]);

% infinite norm on the whole vector
normInf = norm(close_sol - x, inf);
disp('Infinite norm |close_sol - computed_sol|');
disp(normInf);

disp('matrix build time');
disp(matrix_build_time);

disp('solving system time');
disp(solving_system_time);

end