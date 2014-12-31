clear;
clc;

dim = [3, 3];

% build the matrix
[A, b] = build_two_dimensional_problem( dim(1), dim(2),... 
    'interior_2D_2ndOrder_centered_5pt', ...
    'dirichlet_boundary');

%-------- display some debug info ------------------
full(A)
full(b')
% -------- display some debug info ------------------
tstart = tic;

numIter = 1e5;
options.nonneg = true;
[x info restart] = sart(A, b', numIter, [], options);        

telapsed = toc(tstart);
%------------------- display some result -----------------------------
size(x)
reshape(x, dim)
numIter
telapsed