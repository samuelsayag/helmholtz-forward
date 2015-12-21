%USAGE Summary of this function goes here
%   provide a functional example

close all; clear all; clc; format shortE;
%------------------------- parameter of the simulation --------------------
% define the domain of the BVP
a0 = 0;
a1 = 1;
% parameters of the BVP
param.k = 20;
param.m = 1e2;
param.h = (a1-a0)./(param.m - 1);
param.interior = 'std';
param.boundary = 'dirichlet';
param.dirichlet.W = @(params,A,b,i) analytic_sol_1D(param.k, a0-param.h);


% %------------------------- dirichlet --------------------------------------
param.dirichlet.E = @(params,A,b,i) analytic_sol_1D(param.k, a1+param.h);
[ func_scheme, param ] = helmholtz_1D_scheme_factory( param );
[A, b] = build_one_dimensional_problem2(param, func_scheme);
x_dirichlet = A\b;

% %------------------------- sommerfeld -------------------------------------
param.dirichlet = rmfield(param.dirichlet, 'E');
param.boundary = 'sommerfeld_new';
[ func_scheme, param ] = helmholtz_1D_scheme_factory( param );
[A, b] = build_one_dimensional_problem2(param, func_scheme);
x_sommerfeld = A\b;

% analytic solution
d = linspace(a0, a1, param.m);
close_sol = analytic_sol_1D(param.k, d);


figure(1)
subplot(3,1,1);
plot(d, real(close_sol));
title('Helmholtz 1D close solution (Real part)');
axis([0 1 -1 1]);

subplot(3,1,2);
plot(d, real(x_dirichlet));
err_dir = norm(close_sol - x_dirichlet, inf);
tit = sprintf('Dirichlet boundaries, k = %u, NormInf = %1.3e', param.k, err_dir);
title(tit);
axis([0 1 -1 1]);

subplot(3,1,3);
plot(d, real(x_sommerfeld));
err_somm = norm(close_sol - x_sommerfeld, inf);
tit = sprintf('Sommerfeld boundary(East), k = %u, NormInf = %1.3e', param.k, err_somm);
title(tit);
axis([0 1 -1 1]);


% % % show a fraction of the vectors of the solution and the computer one
% % tit = {'close sol' 'computed sol' '(last 10 values)'};
% % disp(tit);
% % disp([close_sol(end-9:end) x(end-9:end)]);
% 
% % infinite norm on the whole vector
% normInf = norm(close_sol - x, inf);
% disp('Infinite norm |close_sol - computed_sol|');
% disp(normInf);
% 
% disp('matrix build time');
% disp(matrix_build_time);
% 
% disp('solving system time');
% disp(solving_system_time);
