close all; clear variables; clc;format short;

n = 100;
dim = [n, n];
% build the matrix
[Af, bf] = build_two_dimensional_problem( dim(1), dim(2),...
    'interior_2D_2ndOrder_centered_5pt', ...
    'dirichlet_boundary');

%-------- display some debug info ------------------
% full(A)
% full(b')
%-------- display some debug info ------------------


tstart = tic;

% solve Ax = b with GAUSS-SEIDEL Method
x = randn(1, size(bf,2));
f = @(x, omega) improved_x(Af, bf, x, omega);
[x, numIter, omega] = gauss_seidel( f, x, 10e6, 10e-6);

telapsed = toc(tstart);

%------------------- display some result -----------------------------
size(x)
solf = transpose(reshape(x, dim));
numIter
telapsed

%-------- display some debug info ------------------
% full(A)
% full(b')
% full(sol)
%-------- display some debug info ------------------

x = linspace(0.125,0.375, n);
y = linspace(0.375, 0.125, n);
[X,Y] = meshgrid( x, y );
axis_scale = [0.125 0.375 0.125 0.375, 0, 110];

mesh(X, Y, solf);
axis(axis_scale)
t4 = sprintf('Heat repartition  R = {(x,y) | 0 < x < 0.5, 0 < y < 0.5}\n... u(0,y) = 0, u(x,0) = 0, u(x,0.5) = 200x, and u(0.5,y) = 200y\n precision:%e, num iter:%e', ...
    1e-6, numIter);
title(t4)