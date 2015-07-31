% test that display a graph of the function k(x,y)

close all; clear variables;

a = 10;
b = 9;
c = 10;
[ pSol, k, beta ] = analytic_sol_3( a, b, c );

param.a = 0;
param.b = pi;
param.c = 0;
param.d = pi;
param.n = 200;
param.h = ( param.b - param.a ) / ( param.n - 1 );
k_dis = discreteWrapper(param.a, param.c, param.h, k);
k_x = @(x) k_dis(x,0);

Discr_dom = linspace(1, param.n, param.n);

plot(Discr_dom, k_x(Discr_dom));
axis([1 param.n -5 25]);
title(sprintf('Plot of k(x) for a = %d; b = %d; c = %d.', a, b, c));
ylabel('k(x)');
xlabel('x');
grid on