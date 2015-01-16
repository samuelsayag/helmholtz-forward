function [ sol ] = analytic_sol_1D( k, x )
%ONE_DIMENSIONAL_ANALYTIC_SOL Summary of this function goes here
%   Detailed explanation goes here
    
    sol = exp(1i * k * x)';
    
end

