% test that display a graph of the function u(x,y)

close all; clear variables;

a = 10; 
b = 9;
c = 10;

[ pSol, k, beta ] = analytic_sol_3( a, b, c );

Dx = [0, pi];
Dy = [0, pi];
x = linspace(Dx(1), Dx(2), 200);
y = linspace(Dy(2), Dy(1), 200);
[X,Y] = meshgrid( x, y );
u = pSol(X, Y);
mesh(X, Y, u);
axis([Dx, Dy, [-0.8 0.8]]);
t1 = sprintf('Plot of u(x) for a = %d; b = %d; c = %d.', a, b, c);
title(t1)
xlabel('x axis');ylabel('y axis'); zlabel('helmholtz');


