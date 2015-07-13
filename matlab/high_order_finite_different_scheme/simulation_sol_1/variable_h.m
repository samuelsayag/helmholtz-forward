% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

% clear all;
close all; clear variables; clc;

% modeled solution
theor = @(x, y, k) helm_sol1_2D( x, y, k );

% basic parameter of the simulation
param.k = 10;
% definition of the place
param.a = 0; 
param.b = 1;
param.c = -1/2; 
param.d = 1/2;

% dirichlet function
param.dirichlet = @(x,y) theor( x, y, param.k );

% define the solver for the Ax=b equation
solver = @(A, b) A\b;

size_g = {12 22 32 47 52 57 62};
error = {1,size(size_g, 2)};
steps = {1,size(size_g, 2)};

for i = 1:size(size_g, 2)    
    param.m = size_g{i};
    param.n = size_g{i};
    param.h = (param.d - param.c)/(param.m-1);
    steps{i} = param.h;
    scheme = Ord6thHelmholtz2D(param.k, param.h);    
    
    ps = ProblemSolver(param, scheme, solver);
    [ A, b, sol ] = ps.solve();
    
    err_cur = ErrorHandler( param, sol );
    error{i} = err_cur.total;
end

r = [size_g; steps; error];
r = [{'size of grid' 'h' 'error'}', r];
mh = cell2mat(steps);
merr = cell2mat(error);
p = polyfit(log(mh), log(merr), 1);

reg = log(mh) * p(1) + p(2);

plot(log(mh), log(merr), log(mh), reg);
title(sprintf('Linear regression of log(err) = log(steps), C=%d', p(2)));
xlabel('log(steps)');
ylabel('log(error)');
legend('log(steps)','log(error)', 'Location', 'southeast'); 

% plot(log(cell2mat(steps)), log(cell2mat(error)));