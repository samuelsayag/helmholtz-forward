close all; clear all; clc;format short;

n = 3;
dim = [n, n];
% build the matrix
[A, b] = build_two_dimensional_problem( dim(1), dim(2),...
    'interior_2D_2ndOrder_centered_5pt', ...
    'dirichlet_boundary');

%-------- display some debug info ------------------
% full(A)
% full(b')
%-------- display some debug info ------------------


tstart = tic;

% solve Ax = b with GAUSS-SEIDEL Method
x = randn(1, size(b,2));
f = @(x, omega) improved_x(A, b, x, omega);
[x, numIter, omega] = gauss_seidel( f, x, 10e6, 10e-6);

telapsed = toc(tstart);

%------------------- display some result -----------------------------
size(x)
sol = transpose(reshape(x, dim));
numIter
telapsed

%-------- display some debug info ------------------
full(A)
full(b')
full(sol)
%-------- display some debug info ------------------


x = linspace(0,0.5, n);
y = linspace(0.5, 0, n);
[X,Y] = meshgrid( x, y );
mesh(X, Y, sol);
axis([0 0.5 0 0.5 0 100])
t4 = sprintf('Heat repartition  R = {(x,y) | 0 < x < 0.5, 0 < y < 0.5}\n... u(0,y) = 0, u(x,0) = 0, u(x,0.5) = 200x, and u(0.5,y) = 200y\n precision:%e, num iter:%e', ...
    1e-6, numIter);
title(t4)