function [A, b, x, n, m] = new_scheme_dirichlet()
%CLASSICAL_SCHEME_DIRICHLET Summary of this function goes here
%   Detailed explanation goes here

A = []; 
b = [];
x = [];

% parameters necessary to compute interior points
k = 10;          
h = 0.5./k;
a = 0; 
b = 1;
n = (b-a)./h;
m = n;

[k,h,a,b,m,n]

% function that will handle the interior points
func_inter = @(A, n, m, i , j, lab) interior_2D_2ndOrder_centered_5pt(...
    A, n, m, i , j, lab, k, h);
% parameters necessary to compute boundary points
S_dirichlet = 1;
W_dirichlet = 0;
N_dirichlet = 0;
E_dirichlet = 0;
% function that will handle the bouindary points
func_bound = @(b, n, m, i, j, lab) dirichlet_boundary(b, n, m, i, j,...
    lab, S_dirichlet, W_dirichlet, N_dirichlet, E_dirichlet);

% create the matrix of finite difference
[A, b] = build_two_dimensional_problem(n, m, ...
    func_inter, func_bound);

%debug
% full(A);
% b;

% solve the system
tstart = tic;
numIter = 1e6;
options.nonneg = true;
[x info restart] = sart(A, b', numIter, [], options);        
telapsed = toc(tstart);

%------------------- display some result -----------------------------
size(x)
reshape(x,n,m)
numIter
telapsed
