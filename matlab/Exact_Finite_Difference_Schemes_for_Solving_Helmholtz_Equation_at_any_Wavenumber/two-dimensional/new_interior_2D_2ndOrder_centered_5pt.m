function [ A ] = new_interior_2D_2ndOrder_centered_5pt(A, n, m, i , j,... 
lab_func, k, h)
% INTERIOR_2D_2NDORDER_CENTERED_5PT 
%  
% return:
%   A: the finite difference matrix
% parameters:
%   A: the finite difference matrix
%   n, m: the born of the system matrix
%   i, j: indice of route in the system matrix
%   lab_func: label creation function for the line in the finite difference
%   matrix.
%   k: [k = omega/c] (omega = 2 * pi * freq, c = speed of sound)
%   h: the step size of our system

l = feval(lab_func, i, j);

% central point
A(l, l) = 4 * besselj(0, k * h);
% south point
if j-1 > 0
    A(l, feval(lab_func, i, j-1)) = -1;
end
% west point
if i-1 > 0
    A(l, feval(lab_func, i-1, j)) = -1;
end
% north point
if j+1 < m+1
    A(l, feval(lab_func, i, j+1)) = -1;
end
% east point
if i+1 < n+1
    A(l, feval(lab_func, i+1, j)) = -1;
end

end

