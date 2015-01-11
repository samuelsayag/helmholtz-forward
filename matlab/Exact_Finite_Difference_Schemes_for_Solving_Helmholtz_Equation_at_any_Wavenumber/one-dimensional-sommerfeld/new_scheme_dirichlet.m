function [A, b, x] = new_scheme_dirichlet()
%CLASSICAL_SCHEME_DIRICHLET Summary of this function goes here
%   Detailed explanation goes here

% parameters necessary to compute interior points
k = 10;
h = 0.5./k;
a = 0;
b = 1;
m = (b-a)./h;
alpha = 1;

[k, h, m, alpha]

func_scheme = @(A, b, i) new_dirichlet_sommerfeld_boundary(A, b, i, k, h,...
    alpha, false);

% create the matrix of finite difference
[A, b] = build_one_dimensional_problem2(m, func_scheme);

% debug
% full(A)
% b

% solve the system
tstart = tic;
numIter = 1e5;
options.nonneg = true;
[x info restart] = sart(A, b', numIter, [], options);        
telapsed = toc(tstart);

%------------------- display some result -----------------------------
size(x)
x
numIter
telapsed
