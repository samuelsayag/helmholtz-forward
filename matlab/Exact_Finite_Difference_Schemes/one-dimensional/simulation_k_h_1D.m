function [ sols, params ] = simulation_k_h_1D( k, h, sim_param )
%SIMULATION_K_H
% perform a simulation for a vector k and h and some fixed parameters
% given by sim_param
% return solution cell vector

dh = size(h,2);
dk = size(k,2);

sols = cell(dh * dk, 1);
params = cell(dh * dk, 1);

cpt_sol = 1;
for i = 1:size(h,2)
    for j = 1:size(k,2)
        sim_param.h = h(1,i);
        sim_param.k = k(1,j);
        sim_param.m = (sim_param.b - sim_param.a)./sim_param.h;
        [ func_scheme, sim_param ] = helmholtz_1D_scheme_factory( sim_param );
        [A, b] = build_one_dimensional_problem2(sim_param, func_scheme);
          sols{cpt_sol, 1} = A\b;
%         sols{cpt_sol, 1} = qmr(A, b, 1e-6, 1e7);          
%         sols{cpt_sol, 1} = bicgstab(A, b, 1e-12, 1e8); % best iterative
%         sols{cpt_sol, 1} = gmres(A,b, 50, 1e-12, 200);
        params{cpt_sol, 1} = sim_param;
        cpt_sol = cpt_sol + 1;
    end
end

end