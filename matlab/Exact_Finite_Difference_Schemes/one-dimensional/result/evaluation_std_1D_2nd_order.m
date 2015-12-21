%USAGE Summary of this function goes here
%   provide a functional example

close all; clear all; clc; format shortE;
%------------------------- parameter of the simulation --------------------
% define the domain of the BVP
a0 = 0;
a1 = 1;
% define gradually larger grid we use
ms = 2.^linspace(3,12,10);
error = zeros(size(ms));
hs = zeros(size(ms));
% parameters of the BVP
param.k = 10;
param.interior = 'std';
% param.boundary = 'sommerfeld_new';
param.boundary = 'dirichlet';


cpt = 1;
for m = ms
    param.m = m;
    param.h = (a1-a0)./(param.m - 1);
    param.dirichlet.W = @(params,A,b,i) analytic_sol_1D(param.k, a0-param.h);
    param.dirichlet.E = @(params,A,b,i) analytic_sol_1D(param.k, a1+param.h);
    
    [ func_scheme, param ] = helmholtz_1D_scheme_factory( param );
    [A, b] = build_one_dimensional_problem2(param, func_scheme);
    
    x = A\b;
    
    d = linspace(a0, a1, param.m);
    close_sol = analytic_sol_1D(param.k, d);
    error(cpt) = norm(close_sol - x, inf);
    hs(cpt) = param.h;
    
    cpt =  cpt + 1;
end

error = error(:);
hs = hs(:);

% fitting of the points
[fitobj, gof, output] = fit(log(hs), log(error),'poly1');
c = coeffvalues(fitobj);

plot(log(hs), log(error), log(hs), log(hs) * c(1) + c(2));
title('Helmholtz 1D BVP: log(error) = log(h)');
legend('log(error) = log(h)', 'linear reg', 'Location', 'northwest');
xlabel('log(h)');
ylabel('log(error)');

fitobj
gof
output