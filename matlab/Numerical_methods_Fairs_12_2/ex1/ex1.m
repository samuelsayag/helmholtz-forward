clear;
clc;

dim = [4, 4];

% build the matrix
[A, b] = build_square_grid( dim(1), dim(2), 'interior_2D_2ndOrder_centered_5pt', ...
    'ex1_boundary');

% debug
% full(A)
% full(b')

% random vector to start with
x = randn(1,size(A,2));

% specific function called to compute line of Gauss-Seidel
f = @(x, omega) improved_x(A, b, x, omega);

tstart = tic;
% sole Ax = b
[x, numIter, omega] = gauss_seidel( f, x, 1000000, 10e-9);

telapsed = toc(tstart);


size(x)
reshape(x, dim - 1)
numIter
telapsed


