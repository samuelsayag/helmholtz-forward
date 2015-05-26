% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

clear variables; close all; clc;

% basic parameter of the simulation
param.k = 10;
param.h = 0.01;
% definition of the area we simulate in it
param.a = 0; 
param.b = 1;
param.c = -1/2; 
param.d = 1/2;
param.m = (param.d - param.c)/param.h + 1;
param.n = (param.b - param.a)/param.h + 1;

% dirichlet function
param.dirichlet = @(x,y) helm_sol1( x, y, param.k );
theta = pi/4; % exact theta calculation
scheme = ExactScheme2D(param.k, param.h, theta);

% define the solver
solver = @(A, b) A\b;

ps = ProblemSolver(param, scheme, solver);
[ A, b, sol ] = ps.solve();

[err, err_r, err_i] = ErrorHandler( param, sol );

param
err
err_r
err_i

% graphical representation
x = linspace(param.a,param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid( x, y );

figure(1)
subplot(1, 2, 1);
mesh(X, Y, real(sol));
title 'Real part'
subplot(1, 2, 2);
mesh(X, Y, imag(sol));
title 'Img part'

figure(2)
theor = helm_sol1( X, Y, param.k );
subplot(1, 2, 1);
mesh(X, Y, real(theor));
title 'Real part'
subplot(1, 2, 2);
mesh(X, Y, imag(theor));
title 'Img part'