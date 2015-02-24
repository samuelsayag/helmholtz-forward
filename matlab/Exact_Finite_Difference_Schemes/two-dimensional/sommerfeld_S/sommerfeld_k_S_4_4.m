%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\..\..\matlab'));
pause on;

% generic parameters of the simulations
sim_param.h = [0.2];
% k = [5];
% k = [150, 100, 70];
% k = sqrt(2)* [30, 25, 20];
k = sqrt(2)* [30, 25, 20, 15, 10, 5];
% k = sqrt(2)* [30, 25, 20];
sim_param.a = 0;
sim_param.b = 0.8;
sim_param.d = 0.8;
sim_param.c = 0;
sim_param.theta = pi/4;
sim_param.m = (sim_param.b - sim_param.a)./sim_param.h;
sim_param.n = (sim_param.d - sim_param.c)./sim_param.h;
% dirichlet boundary
% parameters necessary to compute boundary points
sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j+1) * params.h);
sim_param.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i+1) * params.h, j * params.h);

% declaration of solution structures
sols = {};
params = {};

% simu interior std
sim_param.interior = 'std';
sim_param.boundary = 'std';
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'std';
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior std
sim_param.interior = 'std';
sim_param.boundary = 'new';
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'new';
[ sol, param ] = simulation_k_2D( k, sim_param );
sols = [sols, sol];
params = [params, param];

% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
x = linspace(1, sim_param.m, sim_param.m) * sim_param.h;
y = linspace(sim_param.n, 1, sim_param.n) * sim_param.h;
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
res_kh = cell(size(k,2),1);
res_k = cell(size(k,2),1);
h = sim_param.h;
for i = 1:size(k,2)
    res_kh{i} = h * k(i);
    res_k{i} = k(i);      
end

% % RESULT FOR JUST STANDARD ('std')
% title1 = {'' '' 'E inf' };
% title2 = {'kh' 'k' 'SFD' };
% res_tab = [title1;title2];
% res_tab = [res_tab; res_kh res_k error ];
% res_tab

% % RESULT FOR JUST STANDARD ('new')
% title1 = {'' '' 'E inf' };
% title2 = {'kh' 'k' 'NFD'};
% res_tab = [title1;title2];
% res_tab = [res_tab; res_kh res_k error ];
% res_tab

% RESULT FOR STANDARD AND NEW ('std', 'new')
% just the central scheme
title1 = {'' '' 'SBC' 'SBC' 'NBC' 'NBC'};
title2 = {'kh' 'k' 'SFD' 'NFD' 'SFD' 'NFD' };
res_tab = [title1;title2];
res_tab = [res_tab; res_kh res_k error ];
res_tab


% figure(1)
% x = linspace(1, sim_param.m, sim_param.m) * sim_param.h;
% y = linspace(sim_param.n, 1, sim_param.n) * sim_param.h;
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
%         title 'computed'
%         analytic_tmp = analytic_sol_2D(params{j,i}.k, params{j,i}.theta, X, Y);
%         subplot(2,3, 3+k)
%         plot3(X, Y, real(analytic_tmp));        
%         title 'analytic'
%     end
%     cptFigure = cptFigure + 1;
% end

% pause
% close all;

pause off;