%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\matlab'));

% generic parameters of the simulations
h = [0.02];
k = [5 10 15 20 25 30];
theta = pi/3;
a = 0;
b = 1;
d = 1;
c = 0;
m = (b - a)./h;
n = (d - c)./h;

% time elapsed
tic
% prepare the meshgrid to calculate the analytic solution or to propose
% graphical representation of the solutions
x = linspace(1, m, m) * h;
y = linspace(n, 1, n) * h;
[X,Y] = meshgrid(x,y);

c_k = mat2cell(k', ones(size(k)));
f = @(t) analytic_sol_2D(t, theta, X, Y);
analytic = cellfun(f, c_k,...
    'UniformOutput', false);

% time elapsed
elapsed = toc;
elapsed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICAL REPRESENTATION - BEGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cptFigure = 1;
for i = 1:size(k, 2)
    figure(cptFigure);        
    cpt_grph = mod(i-1, 6) + 1;
    
    subplot(2, 3, cpt_grph)
    plot3(X, Y, real(analytic{i}));            
%     title (sprintf('k: %0.4f',k(i)));
    title (sprintf('k: %u',k(i)));
    
    if cpt_grph == 6
        cptFigure = cptFigure + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICAL REPRESENTATION - END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%