%==========================================================================
% replay the simulation of the article that originated the new schemes for
% 1D and 2D problem.
% These simulations are for the 2D problem
%==========================================================================
%close all; clear all; clc;
pause on;

% generic parameters of the simulations
h = [2e-2];
k = sqrt(2)* [30, 25, 20, 15, 10, 5];
sim_param.a = 0;
sim_param.b = 1;
sim_param.d = 1;
sim_param.c = 0;
sim_param.theta = pi/4;

% declaration of solution structures
sols = {};
params = {};

% simu interior NEW, Sommerfeld boundary NEW
sim_param.interior = 'std';
sim_param.boundary = 'new';
[ sol, param ] = simulation_k_h( k, h, sim_param );
sols = [sols, sol];
params = [params, param];

% simu interior NEW, Sommerfeld boundary NEW
sim_param.interior = 'std';
sim_param.boundary = 'new';
[ sol, param ] = simulation_k_h( k, h, sim_param );
sols = [sols, sol];
params = [params, param];


figure(1)
a = sim_param.a;
b = sim_param.b;
d = sim_param.d;
c = sim_param.c;
h = h(1);
x = linspace(a,b, (b-a)/h + 1);
y = linspace(d,c, (d-c)/h + 1);
[X,Y] = meshgrid(x,y);

cptFigure = 0;
% real part
for i = 1:size(sols,2)
    for j = 1:size(sol,1)
        if j < 4
            figure(2 * cptFigure + 1)
            k = j;
        else
            figure(2 * cptFigure + 2)
            k = j-3;
        end
        subplot(2,3, k)
        plot3(X, Y, real(sols{j,i}));        
        analytic = analytic_sol_2D(params{j,i}.k, params{j,i}.theta, X, Y);
        subplot(2,3, 3+k)
        plot3(X, Y, real(analytic));        
    end
    cptFigure = cptFigure + 1;
end

% pause
% close all;

% error = cell(size(sols));
% analytic = cell(size(sols));
% 
% for i = 1:size(error,1)
%     for j = 1:size(error,2)
%         k = params{i,j}.k;
%         tmp = analytic_sol_2D(k, params.theta, X, Y);
%         analytic{i,j} = tmp(2:end);
%         error{i,j} = norm((analytic{i,j} - sols{i,j}), Inf );
%     end    
% end



% close all;
% 
% mt = sprintf('\t');
% kkk = cell(size(kk));
% for index = 1:length(kk)
%     kkk{index} = kk(index);
% end
% disp (['h = ' num2str(h)])
% % str_vec = {'kh' 'k' 'SFD' 'NFD' 'SFD' 'NFD'};
% str_vec = {'k' 'SFD' 'NFD' 'SFD' 'NFD'};
% error1 = [str_vec; kkk' error]

% tmp = [cell2mat(analytic(1,1)) cell2mat(sols(1,1)) cell2mat(sols(1,4))]

pause off;