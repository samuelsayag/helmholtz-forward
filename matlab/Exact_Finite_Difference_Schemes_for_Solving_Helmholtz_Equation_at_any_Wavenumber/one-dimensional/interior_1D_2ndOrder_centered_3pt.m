function [ A ] = interior_1D_2ndOrder_centered_3pt(A, i, k, h)
% INTERIOR_1D_2NDORDER_CENTERED_3PT 
% Create the adequate line of a point for a 1D second order cented 3 points 
% finite difference method at indice i.
% return: 
%   A: the finiter difference matrix that represent the sytem we determine
% parameter:
%   A: the finite difference matrix that represent the sytem we determine
%   i: the indice of the line we wish to compute
%   k: [k = omega/c] (omega = 2 * pi * freq, c = speed of sound)
%   h: the step size of our system

m = size(A,2);



% A(i, i) = (k * h).^2 - 2;
A(i, i) = 2 - (k * h).^2;

if i-1 > 0
    A(i, i-1) = -1;
end

if i+1 <= m
    A(i, i+1) = -1;
end

end

