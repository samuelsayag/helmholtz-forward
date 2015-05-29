function [ err ] = ErrorHandler( param, computed_sol )
%ERRORHANDLER helper to calculate the error
% param: parameter of the simulation
% computer_sol: the solution that is computed

x = linspace(param.a, param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid( x, y );
analytic = param.dirichlet( X, Y );

diff = abs(transpose(analytic - computed_sol));

err.total = max(sum(diff));
err.real = norm( abs(transpose(real(analytic) - real(computed_sol))));
err.img = norm( abs(transpose(imag(analytic) - imag(computed_sol))));
err.norm = mean(sum(diff));
err.std = std(sum(diff));
err.min = min(sum(diff));

end

