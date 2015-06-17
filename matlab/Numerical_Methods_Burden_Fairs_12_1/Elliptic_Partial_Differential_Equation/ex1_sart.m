close all; clear all; clc;format short;

dim = [50, 50];

% build the matrix
[A, b] = build_two_dimensional_problem( dim(1), dim(2),... 
    'interior_2D_2ndOrder_centered_5pt', ...
    'dirichlet_boundary');

%-------- display some debug info ------------------
% full(A)
% full(b')
% -------- display some debug info ------------------
tstart = tic;

numIter = 1e5;
options.nonneg = true;
[x info restart] = sart(A, b', numIter, [], options);        

telapsed = toc(tstart);
%------------------- display some result -----------------------------
size(x)
sol = transpose(reshape(x, dim));
numIter
telapsed


x = linspace(0,0.5, 50);
y = linspace(0.5, 0, 50);
[X,Y] = meshgrid( x, y );
mesh(X, Y, sol);
axis([0 0.5 0 0.5 0 100])
t4 = sprintf('Heat repartition  R = {(x,y) | 0 < x < 0.5, 0 < y < 0.5}\n u(0,y) = 0, u(x,0) = 0, u(x,0.5) = 200x, and u(0.5,y) = 200y');
title(t4)