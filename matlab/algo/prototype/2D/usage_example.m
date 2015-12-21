function [ A, b, sol, param, close_sol  ] = usage_example()
%USAGE Summary of this function goes here
%   illustration of the usage of the 2D prototype

close all; clear variables; clc;
%------------------------- parameter of the simulation --------------------
% parameter of the equation
param.k = 10;
% define the domain of the BVP
param.h = 1e-2;
param.a = 0;
param.b = 1;
param.d = 1;
param.c = 0;
% compute the grid points
param.m = (param.b - param.a)./param.h +1;
param.n = (param.d - param.c)./param.h +1;
% parameter of the scheme
param.theta = pi/4;
param.bessel = @(x) bessel_exact_theta(x, param.theta);
param.interior = 'new';
param.boundary = 'new';

% definition of the Dirichlet boundaries for the Souht and West side
param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, (j-2) * params.h);
param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-2) * params.h, (j-1) * params.h);

%------------------------- use of the factory -----------------------------
% build the simulation function
[ func_scheme, param ] = helmholtz_2D_scheme_factory( param );

%------------------------- building the matrix ----------------------------
% generate the matrix and vector of the problem
tic 
[A, b] = build_two_dimensional_problem2(param, func_scheme);
matrix_build_time = toc;

%------------------------- solving the system -----------------------------;
sol = A\b;
sol = transpose(reshape(sol, param.m, param.n));
solving_system_time = toc;

%------------------------- checking output --------------------------------
% building the close solution value matrix on the computed domain
x = linspace(param.a, param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid(x,y);
close_sol = analytic_sol_2D(param.k, param.theta, X, Y);

% show for instance extract of the South-East corner values
tit = {'close sol                  ',...
    'computed sol' '(last 2 column 10 values)'};
disp(tit);
cs = full(close_sol); % rename variable for ease reading
s = full(sol); % rename variable for ease reading
disp([cs(size(cs,1)-11:end, size(cs,2)-1:end) ...
       s(size(s,1)-11:end, size(s,2)-1:end)]);

% infinite norm on the whole vector
normInf = norm(cs - s, inf);
disp('Infinite norm of |close_sol - computed_sol|');
disp(normInf);


disp('matrix build time');
disp(matrix_build_time);

disp('solving system time');
disp(solving_system_time);

end


    'close sol           ...'    'computed sol'    '(last 2 column 10 va...'

   0.0757 + 0.9971i   0.0051 + 1.0000i   0.0757 + 0.9971i   0.0051 + 1.0000i
   0.1460 + 0.9893i   0.0757 + 0.9971i   0.1460 + 0.9893i   0.0757 + 0.9971i
   0.2155 + 0.9765i   0.1460 + 0.9893i   0.2155 + 0.9765i   0.1460 + 0.9893i
   0.2840 + 0.9588i   0.2155 + 0.9765i   0.2840 + 0.9588i   0.2155 + 0.9765i
   0.3510 + 0.9364i   0.2840 + 0.9588i   0.3510 + 0.9364i   0.2840 + 0.9588i
   0.4163 + 0.9092i   0.3510 + 0.9364i   0.4163 + 0.9092i   0.3510 + 0.9364i
   0.4795 + 0.8775i   0.4163 + 0.9092i   0.4795 + 0.8775i   0.4163 + 0.9092i
   0.5403 + 0.8415i   0.4795 + 0.8775i   0.5403 + 0.8415i   0.4795 + 0.8775i
   0.5984 + 0.8012i   0.5403 + 0.8415i   0.5984 + 0.8012i   0.5403 + 0.8415i
   0.6535 + 0.7569i   0.5984 + 0.8012i   0.6535 + 0.7569i   0.5984 + 0.8012i
   0.7053 + 0.7089i   0.6535 + 0.7569i   0.7053 + 0.7089i   0.6535 + 0.7569i
   0.7537 + 0.6573i   0.7053 + 0.7089i   0.7537 + 0.6573i   0.7053 + 0.7089i

Infinite norm of |close_sol - computed_sol|
   9.1795e-13

matrix build time
   11.7190

solving system time
   11.7964