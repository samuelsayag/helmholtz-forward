function [ A, b ] = build_one_dimensional_problem( m, func_inter, func_bound )
%BUILD_SQUARE_GRID Summary of this function goes here
%   Detailed explanation goes here

b = zeros(1, m);
A = sparse(m, m);
 
for i = 1:m
    A = feval(func_inter, A, i);
    b = feval(func_bound, b, i);
end

end