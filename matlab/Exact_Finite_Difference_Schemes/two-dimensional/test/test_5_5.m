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

% definition of the Dirichlet boundaries
param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, (j-2) * params.h);
param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-2) * params.h, (j-1) * params.h);
% sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i-1) * params.h, (j) * params.h);
% sim_param.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i) * params.h, (j-1) * params.h);

%------------------------- use of the factory -----------------------------
% build the simulation function
[ func_scheme, param ] = helmholtz_2D_scheme_factory( param );

%------------------------- building the matrix ----------------------------
% generate the matrix and vector of the problem
[A, b] = build_two_dimensional_problem2(param, func_scheme);

%------------------------- solving the system -----------------------------;
sol = A\b;
sol = transpose(reshape(sol, param.m, param.n));

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
disp('Infinite norm');
disp(normInf);

% close all; clear variables; clc;
% addpath(genpath('..\..\..\matlab'));
% pause on;
% 
% % generic parameters of the simulations
% % sim_param.h = 1e-2;
% sim_param.h = 0.2;
% sim_param.k = [5];
% % sim_param.k = [30 * sqrt(2)];
% 
% sim_param.a = 0;
% sim_param.b = 1;
% sim_param.d = 1;
% sim_param.c = 0;
% sim_param.m = (sim_param.b - sim_param.a)./sim_param.h +1;
% sim_param.n = (sim_param.d - sim_param.c)./sim_param.h +1;
% sim_param.theta = pi/4;
% % dirichlet boundary
% % parameters necessary to compute boundary points
% sim_param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i-1) * params.h, (j-2) * params.h);
% sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, (i-2) * params.h, (j-1) * params.h);
% % sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
% %     params.theta, (i-1) * params.h, (j) * params.h);
% % sim_param.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
% %     params.theta, (i) * params.h, (j-1) * params.h);
% 
% sim_param.interior = 'new';
% sim_param.boundary = 'new';
% 
% sim_param
% 
% % build the simulation function
% [ func_scheme, sim_param ] = helmholtz_2D_scheme_factory( sim_param );
% % generate the matrix and vector of the problem
% [A, b] = build_two_dimensional_problem2(sim_param, func_scheme);
% % compute the solution of the equation
% solt = A\b;
% solt = transpose(reshape(solt, sim_param.m, sim_param.n));
% 
% x = linspace(sim_param.a, sim_param.b, sim_param.m);
% y = linspace(sim_param.d, sim_param.c, sim_param.n);
% [X,Y] = meshgrid(x,y);
% 
% analytic = analytic_sol_2D(sim_param.k, sim_param.theta, X, Y);
% error = norm((analytic - solt), Inf );
% 
% error


