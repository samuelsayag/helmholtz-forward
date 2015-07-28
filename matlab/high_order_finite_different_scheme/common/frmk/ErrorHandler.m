function [ err ] = ErrorHandler( param, computed_sol )
%ERRORHANDLER helper to calculate the error
% param: parameter of the simulation
% computer_sol: the solution that is computed

x = linspace(param.a, param.b, param.m);
y = linspace(param.d, param.c, param.n);
[X,Y] = meshgrid( x, y );
analytic = param.dirichlet( X, Y );

diff = analytic - computed_sol;
module_diff = abs(diff);

% norm(X, Inf) = max(sum(abs(X')))
err.normInf = norm(analytic - computed_sol, Inf);
% Attention !!! if A is a matrix norm(A) <> norm(A(:)) 
% norm(A) = max(svd(A))
err.l2err = norm(module_diff(:),2) ./ norm(analytic(:),2);

end

