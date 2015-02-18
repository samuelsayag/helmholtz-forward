function [ sols, params ] = simulation_k_2D( k, sim_param )
%SIMULATION_K_H
% perform a simulation for a vector k and h and some fixed parameters
% given by sim_param
% return solution cell vector

% simultation output structures declaration
dh = 1;
dk = size(k,2);
sols = cell(dh * dk, 1);
params = cell(dh * dk, 1);

cpt_sol = 1;

for j = 1:size(k,2)
    % prepare temporary parameters
    sim_param.k = k(1,j);        
    % build the simulation function
    [ func_scheme, sim_param ] = helmholtz_2D_scheme_factory( sim_param );
    % generate the matrix and vector of the problem
    [A, b] = build_two_dimensional_problem2(sim_param, func_scheme);
    % compute the solution of the equation
      tmp_sol = A\transpose(b);
%     tmp_sol = bicgstab(A, transpose(b), 1e-12, 1e6);
%     tmp_sol = gmres(A, transpose(b), 50, 1e-12, 500);
    tmp_sol = reshape(tmp_sol, sim_param.m, sim_param.n)';
    sols{cpt_sol, 1} = tmp_sol;
    params{cpt_sol, 1} = sim_param;
    cpt_sol = cpt_sol + 1;
end


end

