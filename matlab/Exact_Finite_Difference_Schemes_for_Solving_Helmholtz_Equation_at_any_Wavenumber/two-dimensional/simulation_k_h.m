function [ sols, params ] = simulation_k_h( k, h, sim_param )
%SIMULATION_K_H
% perform a simulation for a vector k and h and some fixed parameters
% given by sim_param
% return solution cell vector

% simultation output structures declaration
dh = size(h,2);
dk = size(k,2);
sols = cell(dh * dk, 1);
params = cell(dh * dk, 1);

% parameters necessary to compute boundary points
sim_param.dirichlet.S = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, i * params.h, (j-1) * params.h);
sim_param.dirichlet.W = @(params, A, b, i, j) analytic_sol_2D(params.k,... 
    params.theta, (i-1) * params.h, j * params.h);

cpt_sol = 1;
for i = 1:size(h,2)
    for j = 1:size(k,2)
        % prepare temporary parameters
        sim_param.h = h(1,i);
        sim_param.k = k(1,j);
        sim_param.m = (sim_param.b - sim_param.a)./sim_param.h;
        sim_param.n = (sim_param.d - sim_param.c)./sim_param.h;        
        % build the simulation function
        [ func_scheme, sim_param ] = helmholtz_2D_scheme_factory( sim_param );
        % generate the matrix and vector of the problem
        [A, b] = build_two_dimensional_problem2(sim_param, func_scheme);
        % compute the solution of the equation
        tmp_sol = gmres(A, b');
        %build the solution of the 2D problem with boundary
        sols{cpt_sol, 1} = get_solution_with_boundaries(tmp_sol, sim_param);
        params{cpt_sol, 1} = sim_param;
        cpt_sol = cpt_sol + 1;
    end
end

end

function [sol_boundary] = get_solution_with_boundaries(sol, params)
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
sol = reshape(sol, params.m, params.n)';
sol_boundary = analytic_sol_2D(k_tmp, params.theta, X, Y);
sol_boundary(1:end-1, 2:end) = sol(:,:);

end