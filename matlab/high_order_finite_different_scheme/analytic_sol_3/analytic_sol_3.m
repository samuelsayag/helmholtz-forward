function [ pSol, k, beta ] = analytic_sol_3( a, b, c )
%BETA_SOL_3 compute 3 different elements:
% 1 - the function pointer that give values of the analytical solution 
% parameterized by a, b, c
% 2 - the beta of the function parameterized by a, b, c
% 3 - the function pointer k(x,y) (variable k parameterized by a, b, c)

if b >= a
    error('k_sol_3:argChck', 'parameters: must satisfy a > b ');
end

if b < 0
    error('k_sol_3:argChck', 'parameters: b >= 0 ');
end

% f = @(x) num2str(x);
% k_str =  strcat(f(a), ' -  ', f(b), ' .* sin( ', f(c) , ' .* x )');

beta = sqrt( a.^2 + b.^2 );

k = @(x,y) a - b .* sin(c .* x);

pSol = @(x,y) exp(- k(x,y)./c ) .* sin(beta .* y); 

end

