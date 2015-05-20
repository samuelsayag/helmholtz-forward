function [ err, err_real, err_img ] = ErrorHandler( param, computed_sol )
%ERRORHANDLER helper to calculate the error
% param: parameter of the simulation
% computer_sol: the solution that is computed

x = linspace(param.a, param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid( x, y );
analytic = param.dirichlet( X, Y );

err = norm( analytic - computed_sol, Inf );
err_real = norm( real(analytic) - real(computed_sol), Inf );
err_img = norm( imag(analytic) - imag(computed_sol), Inf );
end

