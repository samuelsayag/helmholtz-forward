%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\..\..\..\matlab'));

% generic parameters of the simulations
sim_param.h = [0.02];
angle_div = 11;
theta = pi/2 * 1/angle_div * linspace(0, angle_div, angle_div + 1);
d_theta = size(theta, 2);

sim_param.k = 30 * sqrt(2);
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
% simu interior, Sommerfeld North
sim_param.dirichlet.S = S;
sim_param.dirichlet.W = W;
sim_param.dirichlet.E = E;
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];
graph_titles = [graph_titles, sprintf('Sommerfeld on %s side', 'North')];

% simu interior, Sommerfeld East
sim_param.dirichlet.S = S;
sim_param.dirichlet.W = W;
sim_param.dirichlet.N = N;
sim_param.dirichlet = rmfield(sim_param.dirichlet, 'E');
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];
graph_titles = [graph_titles, sprintf('Sommerfeld on %s side', 'East')];

% simu interior, Sommerfeld South
sim_param.dirichlet.E = E;
sim_param.dirichlet.W = W;
sim_param.dirichlet.N = N;
sim_param.dirichlet = rmfield(sim_param.dirichlet, 'S');
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];
graph_titles = [graph_titles, sprintf('Sommerfeld on %s side', 'South')];

% simu interior, Sommerfeld West
sim_param.dirichlet.E = E;
sim_param.dirichlet.S = S;
sim_param.dirichlet.N = N;
sim_param.dirichlet = rmfield(sim_param.dirichlet, 'W');
[ sol, param ] = simulation_theta_2D( theta, sim_param );
sols = [sols, sol];
params = [params, param];
graph_titles = [graph_titles, sprintf('Sommerfeld on %s side', 'West')];

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

% % FOUR COLUMN RESULT
% just the central scheme
title1 = {'theta' 'N' 'E' 'S' 'W' };
res_cos_theta = mat2cell(cos(theta'), ones(1, d_theta));
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
   subplot(2,2, i);
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