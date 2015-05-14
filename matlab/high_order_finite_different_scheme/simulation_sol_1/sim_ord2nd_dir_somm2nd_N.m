% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

clear variables;
close all; clc;

% basic parameter of the simulation
param.k = 20;
param.h = 0.001;
% definition of the area we simulate in it
param.a = 0; 
param.b = 1;
param.c = -1/2; 
param.d = 1/2;
param.m = (param.d - param.c)/param.h + 1;
param.n = (param.b - param.a)/param.h + 1;

% dirichlet function
param.dirichlet = @(x,y) helm_sol1( x, y, param.k );
param.north = 'sommerfeld';
scheme = Ord2ndHelmholtz2D(param.k, param.h);
sommerfeld = Ord2ndSommerfeld2D( param.h, param.k, scheme );
solver = @(A, b) A\b;

ps = ProblemSolver(param, scheme, solver, sommerfeld);
tic
[ A, b, sol ] = ps.solve();
toc            
% full(A)
% full(b)
% full(x)            

[err, err_r, err_i] = ErrorHandler( param, sol );
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