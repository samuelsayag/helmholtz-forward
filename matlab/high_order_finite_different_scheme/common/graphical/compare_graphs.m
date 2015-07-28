function compare_graphs( param,  sol, error, axis_scale, f_num)
%COMPARE_GRAPHS 
% param: a struct that describe the parameters for a simulation as describe
% in the classes: ProblemSolver and BasicScheme
% sol: the computed solution (complex matrix) that we wish to display
% against the theoretical solution.
% error: a struct that has three fields (error.total, error.real,
% error.img).

if (nargin == 3)
    f_num = 1;
end

% graphical representation
x = linspace(param.a,param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid( x, y );

figure(f_num)
subplot(2, 2, 1);
mesh(X, Y, real(sol));
axis(axis_scale)
t1 = sprintf('Computed Solution (Real Part) \nErr: Norm Inf %e, \nrelative %e',...
    error.normInf, error.l2err);
title(t1)
xlabel('x axis');ylabel('y axis'); zlabel('helmholtz');

subplot(2, 2, 2);
mesh(X, Y, imag(sol));
axis(axis_scale)
t2 = sprintf('Computed solution (Imaginary Part) \nErr: Norm Inf %e, \nrelative %e',...
    error.normInf, error.l2err);
title(t2)
xlabel('x axis');ylabel('y axis'); zlabel('helmholtz');

theor = param.dirichlet(X, Y);
subplot(2, 2, 3);
mesh(X, Y, real(theor));
axis(axis_scale)
t3 = sprintf('Closed Solution (Real Part)');
title(t3)
xlabel('x axis');ylabel('y axis'); zlabel('helmholtz');

subplot(2, 2, 4);
mesh(X, Y, imag(theor));
axis(axis_scale)
t4 = sprintf('Closed Solution (Imaginary Part)');
title(t4)
xlabel('x axis');ylabel('y axis'); zlabel('helmholtz');

end

