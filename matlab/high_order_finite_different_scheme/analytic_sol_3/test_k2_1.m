% test that display a graph of the function k(x,y)

close all; clear variables;

a = 10; 
b = 9;
c = 10;

[ pSol, k, beta ] = analytic_sol_3( a, b, c );
k2 = @(x,y) k(x,y).^2;

Dx = [0, pi];
x = linspace(Dx(1), Dx(2), 200);
y = 1; % do not use it but we want k as a general function of x and y

plot(x, k2(x,y));
axis([Dx, [-5 400]])
title(sprintf('Plot of k(x) for a = %d; b = %d; c = %d.', a, b, c));
ylabel('k(x)');
xlabel('x');
grid on