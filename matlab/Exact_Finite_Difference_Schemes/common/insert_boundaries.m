 function [ sol_boundary ] = insert_boundaries( sol, params )
%INSERT_BOUNDARIES Summary of this function goes here
%   add boundaries to a calculated solution
% the analytic solution is precalculated and will server two purpose:
% 1 - calculation of the error norm((A_analytic - A_computed),inf)
% 2 - the middle value are replaced by the computed one to build the whole
% solution so that the dirichlet may be kept.
% instead of computing just the dirichlet value  
x = linspace(params.a,params.b, params.m+1);
y = linspace(params.d,params.c, params.n+1);
[X,Y] = meshgrid(x,y);
% build experimental sol WITH ITS BOUNDARY
k_tmp = params.k;
sol_boundary = analytic_sol_2D(k_tmp, params.theta, X, Y);
sol_boundary(1:end-1, 2:end) = sol(:,:);

end

