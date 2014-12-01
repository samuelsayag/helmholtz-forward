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

% solve Ax = b with GAUSS-SEIDEL Method
x = randn(1, size(b,2));
f = @(x, omega) improved_x(A, b, x, omega);
[x, numIter, omega] = gauss_seidel( f, x, 10e6, 10e-6);

telapsed = toc(tstart);

%------------------- display some result -----------------------------
size(x)
reshape(x, dim)
numIter
telapsed


