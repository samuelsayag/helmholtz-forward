function [A, b, x] = classical_scheme_dirichlet()
%CLASSICAL_SCHEME_DIRICHLET Summary of this function goes here
%   Detailed explanation goes here

% parameters necessary to compute interior points
k = 10;
h = 0.5./k;
a = 0;
b = 1;
m = (b-a)./h;
func_inter = @(A, i) interior_1D_2ndOrder_centered_3pt(A, i, k, h);

% parameters necessary to compute boundary points
alpha = 0;
beta = 0.1;
func_bound = @(b, i) dirichlet_boundary(b, i, alpha, beta);

% create the matrix of finite difference
[A, b] = build_one_dimensional_problem(m, func_inter, func_bound);

debug
full(A)
b

% % solve the system
% tstart = tic;
% numIter = 1e6;
% options.nonneg = true;
% [x info restart] = sart(A, b', numIter, [], options);        
% telapsed = toc(tstart);
% 
% %------------------- display some result -----------------------------
% size(x)
% x
% numIter
% telapsed
