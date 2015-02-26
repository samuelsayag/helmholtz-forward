%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
close all; clear all; clc;
addpath(genpath('..\..\..\..\matlab'));

% generic parameters of the simulations
h = [0.02];
angle_div = 6;
theta = pi/2 * 1/angle_div * linspace(0, angle_div, angle_div + 1);
k = 15 * sqrt(2);
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

dimTheta = size(theta, 2);
analytic = cell(size(dimTheta, 2));
for i = 1:dimTheta
    analytic{i} = analytic_sol_2D(k, theta(i), X, Y);    
end

% time elapsed
elapsed = toc;
elapsed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICAL REPRESENTATION - BEGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cptFigure = 1;

for i = 1:dimTheta
    figure(cptFigure);        
    k = mod(i-1, 6) + 1;
    subplot(2, 3, k)
    plot3(X, Y, real(analytic{i}));        
    title (sprintf('theta: %0.4f',theta(i)));
    if k == 6
        cptFigure = cptFigure + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICAL REPRESENTATION - END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%