clear;
clc;

dim = [3, 3];

% build the matrix
[A, b] = build_rect_grid( dim(1), dim(2), 'interior_2D_2ndOrder_centered_5pt', ...
    'ex1_boundary');

%-------- display some debug info ------------------
% full(A)
% full(b')
%-------- display some debug info ------------------
tstart = tic;

numIter = 200;
[x info restart] = sart(A, b', numIter);        

telapsed = toc(tstart);
%------------------- display some result -----------------------------
size(x)
reshape(x, dim)
numIter
telapsed