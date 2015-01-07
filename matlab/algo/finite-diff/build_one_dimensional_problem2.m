function [ A, b ] = build_one_dimensional_problem2( m, func_scheme)
%BUILD_SQUARE_GRID Summary of this function goes here
%   Detailed explanation goes here

b = zeros(1, m);
A = sparse(m, m);
 
for i = 1:m
    [ A , b ] = feval(func_scheme, A, b, i);
end

end