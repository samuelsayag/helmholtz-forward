%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\..\..\matlab'));
pause on;

% generic parameters of the simulations
h = [2e-2];
k = sqrt(2)* [30, 25, 20, 15, 10, 5];
sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.theta = pi/4;

% declaration of solution structures
sols = {};
params = {};

% simu interior NEW, Sommerfeld boundary NEW
sim_param.interior = 'std';
sim_param.boundary = 'std';
[ sol, param ] = simulation_k_h_2D( k, h, sim_param );
sols = [sols, sol];
params = [params, param];

% % simu interior NEW, Sommerfeld boundary NEW
% sim_param.interior = 'std';
% sim_param.boundary = 'new';
% [ sol, param ] = simulation_k_h_2D( k, h, sim_param );
% sols = [sols, sol];
% params = [params, param];


% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
a = sim_param.a;
b = sim_param.b;
d = sim_param.d;
c = sim_param.c;
h = h(1);
x = linspace(a,b, (b-a)/h + 1);
y = linspace(d,c, (d-c)/h + 1);
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
    res_kh{i} = h(1) * k(i);
    res_k{i} = k(i);
    exact_theta{i}  = bessel_exact_theta(res_kh{i}, sim_param.theta);
    J0_kh{i} = besselj(0, res_kh{i});       
end

% RESULT FOR JUST STANDARD ('std,'std')
title1 = {'' '' 'E inf'  'J0(kh)' 'J0(kh)'};
title2 = {'kh' 'k' 'SFD' '[0,pi]' 'Exact Theta'};
res_tab = [title1;title2];
res_tab = [res_tab; res_kh res_k error J0_kh exact_theta];
res_tab

% title1 = {'' '' 'E inf' 'E inf' 'J0(kh)' 'J0(kh)'};
% title2 = {'kh' 'k' 'SFD' 'NFD' '[0,pi]' 'Exact Theta'};
% res_tab = [title1;title2];
% res_tab = [res_tab; res_kh res_k error J0_kh exact_theta];
% res_tab


figure(1)
a = sim_param.a;
b = sim_param.b;
d = sim_param.d;
c = sim_param.c;
h = h(1);
x = linspace(a,b, (b-a)/h + 1);
y = linspace(d,c, (d-c)/h + 1);
[X,Y] = meshgrid(x,y);

cptFigure = 0;
% real part
for i = 1:size(sols,2)
    for j = 1:size(sol,1)
        if j < 4
            figure(2 * cptFigure + 1)
            k = j;
        else
            figure(2 * cptFigure + 2)
            k = j-3;
        end
        subplot(2,3, k)
        plot3(X, Y, real(sols{j,i}));        
        analytic = analytic_sol_2D(params{j,i}.k, params{j,i}.theta, X, Y);
        subplot(2,3, 3+k)
        plot3(X, Y, real(analytic));        
    end
    cptFigure = cptFigure + 1;
end

% pause
% close all;

pause off;