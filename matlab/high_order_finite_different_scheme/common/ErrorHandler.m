function [ err, err_real, err_img ] = ErrorHandler( param, computer_sol )
%ERRORHANDLER helper to calculate the error
% param: parameter of the simulation
% computer_sol: the solution that is computed

x = linspace(param.a, param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid( x, y );
analytic = param.dirichlet( X, Y );

err = norm( analytic - computer_sol, Inf );
err_real = norm( real(analytic) - real(computer_sol), Inf );
err_img = norm( imag(analytic) - imag(computer_sol), Inf );
end

