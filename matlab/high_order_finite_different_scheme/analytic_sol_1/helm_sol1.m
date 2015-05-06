function [ helm ] = helm_sol1( X, Y, k )
%HELM_SOL_1 Summary of this function goes here
%   Detailed explanation goes here
    beta = sqrt(k.^2 - pi.^2);
    helm = cos(pi * Y) * exp(-1i * beta * X);     
end

