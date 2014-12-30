function [ A ] = interior_2D_2ndOrder_centered_5pt(A, n, m, i , j, lab_func )
% Create the adequate line of a point for a 2D second order cented 5 points 
% finite difference method.
%INTERIOR_2D_2NDORDER_CENTERED_5PT Summary of this function goes here
%   Detailed explanation goes here

% l = feval(lab_func, i, j, n, m);
l = feval(lab_func, i, j);

A(l, l) = 4;

if i-1 > 0
    A(l, feval(lab_func, i-1, j)) = -1;
end

if j+1 < m
    A(l, feval(lab_func, i, j+1)) = -1;
end

if i+1 < n
    A(l, feval(lab_func, i+1, j)) = -1;
end

if j-1 > 0
    A(l, feval(lab_func, i, j-1)) = -1;
end

end

