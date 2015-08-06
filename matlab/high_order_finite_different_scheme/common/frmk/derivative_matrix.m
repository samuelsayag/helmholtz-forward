function [ f, fx, fy, fxx, fyy ] = derivative_matrix( f_mat, h )
% DERIVATIVE_MATRIX compute the first and second derivative of the matrix
% with respect to x and to y.
% Due to the type of difference used for the derivative (central difference
% for the first derivative (y, y) and also for the second derivative, the 
% matrix given must be at least [3x3].
% If the original matrix was [m x n], the computed result is list of matrix
% that size are [m-1 x n-1].
% The original matrix is also returned with its new reduced dimensions.
% The error due to the mathematical technic is in O(h.^2) for all the 
% derivatives.
%
% f_mat: the matrix of the function that will be derived ([m x n])
% h : the length of the step between the points, it must be equal along x
% and y;

f = f_mat(2:end-1, 2:end-1);

% -------------------------------------------------------------------------
% computation of the elementary differences
% -------------------------------------------------------------------------

fxplush = f_mat(2:end-1, 3:end);
fxminush = f_mat(2:end-1, 1:end-2);

fyplush = f_mat(1:end-2, 2:end-1);
fyminush = f_mat(3:end, 2:end-1);

% -------------------------------------------------------------------------
% computation of the first derivative with respect to X
% -------------------------------------------------------------------------
% central difference O(h^2)
fx = (fxplush - fxminush)  ./ (2*h);

% -------------------------------------------------------------------------
% computation of the first derivative with respect to X
% -------------------------------------------------------------------------
% central difference O(h^2)
fy = (fyplush - fyminush)  ./ (2*h);

% -------------------------------------------------------------------------
% computation of the second derivative with respect to X
% -------------------------------------------------------------------------
% computation of the second derivative O(h^2)
fxx = (fxplush - 2 * f + fxminush) ./ (h.^2);

% -------------------------------------------------------------------------
% computation of the second derivative with respect to Y
% -------------------------------------------------------------------------
% computation of the second derivative O(h^2)
fyy = (fyplush - 2 * fy + fyminush) ./ (h.^2);

end

