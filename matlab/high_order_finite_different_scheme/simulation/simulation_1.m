% Problem taken 
% ITERATIVE SCHEMES FOR HIGH ORDER COMPACT DISCRETIZATIONS
% TO THE EXTERIOR HELMHOLTZ EQUATION

clear;clc;

% basic parameter of the simulation
param.k = 10;
param.h = 0.1;

% definition of the place
param.a = 0; 
param.b = 1;
param.c = -1/2; 
param.d = 1/2;
param.m = (param.d - param.c)/param.h + 1;
param.n = (param.b - param.a)/param.h + 1;

% dirichlet function
x_val = @(x) param.a + (x-1) * param.h;
y_val = @(x) param.c + (x-1) * param.h;
param.dirichlet = @(i,j) helm_sol1( x_val(j), y_val(i), param.k );

scheme = Ord2ndHelmholtz2D(param.k, param.h);
solver = @(A, b) A\b;
ps = ProblemSolver(param, scheme, solver);
[ A, b, x ] = ps.solve();
            
full(A)
full(b)
full(x)            

x_sol = linspace(param.a,param.b, param.m);
y_sol = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid(x_sol,y_sol);

figure(1)
subplot(1, 2, 1);
mesh(X, Y, real(x));
title 'Real part'
subplot(1, 2, 2);
mesh(X, Y, imag(x));
title 'Img part'