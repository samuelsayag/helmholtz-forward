%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('matlab'));
pause on;

% generic parameters of the simulations
sim_param.h = [0.02];
angle_div = 11;
theta = pi/2 * 1/angle_div * linspace(0, angle_div, angle_div + 1);
% theta =   1.1610 + abs(1.2976 - 1.1610) * 1/angle_div * linspace(0, angle_div, angle_div + 1);
% theta =   1.2263 + abs(1.2382 - 1.2263) * 1/angle_div * linspace(0, angle_div, angle_div + 1);
d_theta =size(theta,2);

sim_param.k = 30 * sqrt(2);
sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.m = (sim_param.b - sim_param.a)./sim_param.h;
sim_param.n = (sim_param.d - sim_param.c)./sim_param.h;
% dirichlet boundary
% parameters necessary to compute boundary points
% sim_param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
%     params.theta, i * params.h, (j-1) * params.h);
sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
sim_param.dirichlet.N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j+1) * params.h);
sim_param.dirichlet.E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i+1) * params.h, j * params.h);

% declaration of solution structures
sols = {};
params = {};
tic

% simu interior NEW
sim_param.interior = 'std';
sim_param.boundary = 'std';
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'std';
sim_param.boundary = 'new';
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'std';
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW
sim_param.interior = 'new';
sim_param.boundary = 'new';
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];

% time elapsed
elapsed = toc;
elapsed

% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
res_theta = mat2cell(theta', ones(1, d_theta));
x = linspace(1, sim_param.m, sim_param.m) * sim_param.h;
y = linspace(sim_param.n, 1, sim_param.n) * sim_param.h;
[X,Y] = meshgrid(x,y);

% calculate the error for each simulation
error = cell(size(sols));
analytic = cell(size(sols));
fa = @(t) analytic_sol_2D(sim_param.k, t, X, Y);
fe = @(a, s) norm((a - s), Inf );
for j = 1:size(error,2)
    analytic(:,j) = cellfun(fa, res_theta, 'UniformOutput', false);
    error(:,j) = cellfun( fe, analytic(:,j), sols(:,j), 'UniformOutput', false );
end

% % ONE COLUMN RESULT
% % just the central scheme
title1 = {'' '' 'error' 'error' 'error' 'error'};
title2 = {'theta' 'cos(theta)' 'SFD - SBC' 'SFD - NBC' 'NFD - SBC' 'NFD - NBC'};
res_tab = [title1; title2];
res_cos_theta = mat2cell(cos(theta'), ones(1, d_theta));
res_tab = [res_tab; res_theta res_cos_theta error ];
res_tab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN - GRAPHICAL REPRESENTATION  OF THE CURVES error = f(cos theta) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cptFigure = 1;

% preparation of the results
cptFigure = cptFigure + 1;
figure(cptFigure);

nb_col = 2;
for i = 1:size(error,2)
    cpt_graph = mod(i-1, nb_col * 2) + 1;
    subplot(2,2,cpt_graph);    
    plot(cos(theta), cell2mat(error(:,i)));
    title (sprintf('%s', title2{2+i}));
    xlabel 'cos(theta)';
    ylabel 'error';    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END - GRAPHICAL REPRESENTATION  OF THE CURVES error = f(cos theta) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN - GRAPHICAL REPRESENTATION  OF THE CURVES real(Z) = f(X,Y) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(cptFigure)
% x = linspace(1, sim_param.m, sim_param.m) * sim_param.h;
% y = linspace(sim_param.n, 1, sim_param.n) * sim_param.h;
% [X,Y] = meshgrid(x,y);
% 
% nb_col = 3;
% cptFigure = 1;
% for i = 1:d_theta
%     figure(cptFigure); % index of figure    
%     k = mod(i-1, nb_col) + 1; % index of graph in the figure
%     
%     subplot(2, nb_col, k)
%     plot3(X, Y, real(analytic{i}));        
%     title (sprintf('analytic theta: %0.4f',theta(i)));
% 
%     subplot(2, nb_col, nb_col + k)
%     plot3(X, Y, real(sols{i}));        
%     title (sprintf('computed theta: %0.4f, error: %0.4f',theta(i), error{i}));    
%     
%     if k == nb_col && i < d_theta
%         cptFigure = cptFigure + 1;
%     end
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END - GRAPHICAL REPRESENTATION  OF THE CURVES real(Z) = f(X,Y) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%