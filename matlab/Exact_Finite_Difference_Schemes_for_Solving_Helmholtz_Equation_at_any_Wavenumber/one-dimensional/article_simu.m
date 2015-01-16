%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulation are for the 1D problem
%==========================================================================
close all; clear all; clc;

h = [1e-2];
k = [1e1, 3e1, 5e1, 7e1, 1e2, 1.5e2];
sim_param.a = 0;
sim_param.b = 1;
sim_param.dirichlet.W = @(params, A, b, i) 1;

% simu SBC - SFD
sim_param.interior = 'std';
sim_param.boundary = 'sommerfeld_std';
[ sol, param ] = simulation_k_h( k, h, sim_param );
sols = sol;
params = param;

% simu SBC - NFD
sim_param.interior = 'new';
sim_param.boundary = 'sommerfeld_new';
[ sol, param ] = simulation_k_h( k, h, sim_param );
sols = [sols,sol];
params = [params, param];

% simu NBC - SFD
sim_param.interior = 'new';
sim_param.boundary = 'sommerfeld_std';
[ sol, param ] = simulation_k_h( k, h, sim_param );
sols = [sols,sol];
params = [params, param];

% simu NBC - NFD
sim_param.interior = 'new';
sim_param.boundary = 'sommerfeld_new';
[ sol, param ] = simulation_k_h( k, h, sim_param );
sols = [sols,sol];
params = [params, param];

error = cell(size(sols));
analytic = cell(size(sols));

for i = 1:size(error,1)
    for j = 1:size(error,2)
        k = params{i,j}.k;
        x = linspace(0,1, params{i,j}.m );
        analytic{i,j} = analytic_sol_1D(k, x);
        error{i,j} = norm((analytic{i,j} - sols{i,j}), Inf );
    end    
end

i = 2;
figure % create new figure
subplot(2,3,1) % first subplot
plot(analytic{i,1});
title('Analytic sol')
axis equal

subplot(2,3,2) % first subplot
plot(sols{i,1});
title('SBC-SFD')
axis equal

subplot(2,3,3) % first subplot
plot(sols{i,2});
title('SBC-NFD')
axis equal

subplot(2,3,4) % first subplot
plot(sols{i,3});
title('NBC-SFD')
axis equal

subplot(2,3,5) % first subplot
plot(sols{i,4});
title('NBC-NFD')
axis equal


