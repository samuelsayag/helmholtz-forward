%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear variables; clc;
addpath(genpath('..\..\..\matlab'));
pause on;

% generic parameters of the simulations
sim_param.h = 2e-2;
k = sqrt(2)* [30, 25, 20, 15, 10, 5];
% k = sqrt(2)* [5 ];
sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.m = (sim_param.b - sim_param.a)./sim_param.h +1;
sim_param.n = (sim_param.d - sim_param.c)./sim_param.h +1;
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

% simu interior NEW, Sommerfeld boundary NEW
sim_param.interior = 'std';
sim_param.boundary = 'new';
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW, Sommerfeld boundary NEW
sim_param.interior = 'new';
sim_param.boundary = 'new';
sim_param.bessel = @(x) bessel_exact_theta(x, sim_param.theta);
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
x = linspace(sim_param.a, sim_param.b, sim_param.m);
y = linspace(sim_param.d, sim_param.c, sim_param.n);
[X,Y] = meshgrid(x,y);

% calculate the error for each simulation
error = cell(size(sols));
analytic = cell(size(sols));
for i = 1:size(error,1)
    for j = 1:size(error,2)
        k_temp = params{i,j}.k;
        analytic{i,j} = analytic_sol_2D(k_temp, sim_param.theta, X, Y);
        error{i,j} = norm((analytic{i,j} - sols{i,j}), Inf );
    end    
end


% preparation of the results
J0_kh = cell(size(k,2),1);
exact_theta = cell(size(k,2),1);
res_kh = cell(size(k,2),1);
res_k = cell(size(k,2),1);
for i = 1:size(k,2)
    res_kh{i} = sim_param.h(1) * k(i);
    res_k{i} = k(i);
    exact_theta{i}  = bessel_exact_theta(res_kh{i}, sim_param.theta);
    J0_kh{i} = besselj(0, res_kh{i});       
end

% % RESULT FOR JUST STANDARD ('std,'std')
% title1 = {'' '' 'E inf'  'J0(kh)' 'J0(kh)'};
% title2 = {'kh' 'k' 'SFD' '[0,pi]' 'Exact Theta'};
% res_tab = [title1;title2];
% res_tab = [res_tab; res_kh res_k error J0_kh exact_theta];
% res_tab

title1 = {'' '' 'E inf' 'E inf' 'J0(kh)' 'J0(kh)'};
title2 = {'kh' 'k' 'SFD' 'NFD' '[0,pi]' 'Exact Theta'};
res_tab = [title1;title2];
res_tab = [res_tab; res_kh res_k error J0_kh exact_theta];
res_tab


% figure(1)
% x = linspace(params.a, params.b, params.m);
% y = linspace(params.d, params.c, params.n);
% [X,Y] = meshgrid(x,y);
% 
% cptFigure = 0;
% % real part
% for i = 1:size(sols,2)
%     for j = 1:size(sol,1)
%         if j < 4
%             figure(2 * cptFigure + 1)
%             k = j;
%         else
%             figure(2 * cptFigure + 2)
%             k = j-3;
%         end
%         subplot(2,3, k)
%         plot3(X, Y, real(sols{j,i}));        
%         analytic = analytic_sol_2D(params{j,i}.k, params{j,i}.theta, X, Y);
%         axis_scale = [sim_param.a, sim_param.b, sim_param.c, sim_param.d, -1, 1];
%         axis(axis_scale)
%         subplot(2,3, 3+k)
%         plot3(X, Y, real(analytic));        
%         axis_scale = [sim_param.a, sim_param.b, sim_param.c, sim_param.d, -1, 1];
%         axis(axis_scale)
%     end
%     cptFigure = cptFigure + 1;
% end

% pause
% close all;

pause off;