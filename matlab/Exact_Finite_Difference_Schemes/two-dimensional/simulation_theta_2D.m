function [ sols, params, As, Bs ] = simulation_theta_2D( theta, sim_param )
%SIMULATION_K_H
% perform a simulation for a vector k and h and some fixed parameters
% given by sim_param
% return solution cell vector

% simultation output structures declaration
dtheta = size(theta,2);

sols = cell(dtheta, 1);
params = cell(dtheta, 1);

As = cell(dtheta, 1);
Bs = cell(dtheta, 1);

cpt_sol = 1;
sim_param
for j = 1:dtheta
    % prepare temporary parameters
    sim_param.theta = theta(1,j);        
    % build the simulation function
    [ func_scheme, sim_param ] = helmholtz_2D_scheme_factory( sim_param );
    % generate the matrix and vector of the problem
    [A, b] = build_two_dimensional_problem2(sim_param, func_scheme);
    As{cpt_sol, 1} = A;
    Bs{cpt_sol, 1} = b;
    % compute the solution of the equation
    tmp_sol = A\b;    
%     tmp_sol = bicgstab(A, b, 1e-6, 1e5);
%       tmp_sol = qmr(A, b, 1e-6, 1e5);
%     tmp_sol = gmres(A, b, 60, 1e-6, 300);
    tmp_sol = transpose(reshape(tmp_sol, sim_param.m, sim_param.n));
    sols{cpt_sol, 1} = tmp_sol;
    params{cpt_sol, 1} = sim_param;
    cpt_sol = cpt_sol + 1;
end


end

