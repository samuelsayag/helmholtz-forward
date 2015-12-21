function [ g ] = define_grid(r_width, r_height, accuracy, k2 )
% function [ g ] = make_grid(width, height, accuracy)
% we build here a matrix of width and height given as parameter and a
% precision given by accuracy.
% r_width : 
% r_height : 
% accuracy : scalar 
% k2 : scalar

% centered the width and height
r_width = r_width - r_width(1);
r_height = r_height - r_height(1);
% create the step with given precision
h = accuracy / k2;
% create the sparse matrix
m =ceil(r_height(2) / h)
n = ceil(r_width(2) / h)
g = sparse(m, n);

end

