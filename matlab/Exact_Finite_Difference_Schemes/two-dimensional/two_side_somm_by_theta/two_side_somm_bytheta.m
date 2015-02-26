%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\..\..\..\matlab'));

% generic parameters of the simulations
sim_param.h = [0.02];
angle_div = 119;
theta = pi/2 * 1/angle_div * linspace(0, angle_div, angle_div + 1);

sim_param.k = 15 * sqrt(2);
sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.m = (sim_param.b - sim_param.a)./sim_param.h;
sim_param.n = (sim_param.d - sim_param.c)./sim_param.h;
sim_param.interior = 'new';
sim_param.boundary = 'new';

% dirichlet boundary
% parameters necessary to compute boundary points
S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j-1) * params.h);
W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);
N = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j+1) * params.h);
E = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i+1) * params.h, j * params.h);

% declaration of solution structures
sols = {};
params = {};
% just for graphique
graph_titles = {};

tic
% simu interior, Sommerfeld East-North
sim_param.dirichlet.E = E;
sim_param.dirichlet.N = N;
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];
graph_titles = [graph_titles, sprintf('Error = f(cos(theta)),  Sommerfeld on %s side', 'N-E')];

% simu interior, Sommerfeld East-North
sim_param.dirichlet.S = S;
sim_param.dirichlet.W = W;
sim_param.dirichlet = rmfield(sim_param.dirichlet, 'E');
sim_param.dirichlet = rmfield(sim_param.dirichlet, 'N');
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];
graph_titles = [graph_titles, sprintf('Error = f(cos(theta)),  Sommerfeld on %s side', 'S-W')];

% time elapsed
elapsed = toc;
elapsed

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
        theta_temp = params{i,j}.theta;
        analytic{i,j} = analytic_sol_2D(sim_param.k, theta_temp, X, Y);
        error{i,j} = norm((analytic{i,j} - sols{i,j}), Inf );
    end    
end

% preparation of the results
res_theta = cell(size(theta,2),1);
h = sim_param.h;
for i = 1:size(theta,2)
    res_theta{i} = theta(i);      
end

% % FOUR COLUMN RESULT
% just the central scheme
title1 = {'theta' 'N-E' 'S-W' };
res_tab = [title1; res_theta error ];
res_tab


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GRAPHICAL REPRESENTATION - ERROR = f(cos(theta))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

for i = 1:size(sols,2)
errorY = zeros(size(theta,2),1);
cosTheta = zeros(size(theta,2),1);
    for j = 1:size(sols,1)
        errorY(j) = error{j, i};      
        cosTheta(j)=cos(theta(j));        
    end
   subplot(1,2, i);
   plot(cosTheta, errorY);
   title(graph_titles{i});
   xlabel 'cos(theta)'
   ylabel 'Error'
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GRAPHICAL REPRESENTATION - BEGIN
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cptFigure = 1;
% 
% for i = 1:dimTheta
%     figure(cptFigure);        
%     k = mod(i-1, 6) + 1;
%     subplot(2, 3, k)
%     plot3(X, Y, real(analytic{i}));        
%     title (sprintf('theta: %0.4f',theta(i)));
%     if k == 6
%         cptFigure = cptFigure + 1;
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GRAPHICAL REPRESENTATION - END
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%