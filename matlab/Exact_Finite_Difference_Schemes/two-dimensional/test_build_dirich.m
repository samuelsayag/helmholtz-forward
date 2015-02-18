close all; clear all; clc;

param.m = 4;
param.n = 4;
param.k = 30 * sqrt(2);
param.h = 0.25;
param.theta = pi/4;

f = @(params, A, b, i, j) analytic_sol_2D(params.k, params.theta, ...
    i * params.h, j * params.h);

dirichx = zeros(1,param.m);
j = 0;
for i=1:param.m
    dirichx(i) = f(param, [], [], i, j);
end

dirichy = zeros(1,param.n);
i = 0;
for j=1:param.m
    dirichy(j) = f(param, [], [], i, j);
end