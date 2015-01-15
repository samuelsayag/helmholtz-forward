function [ A, b ] = build_one_dimensional_problem2( params, func_scheme)
%BUILD_SQUARE_GRID Summary of this function goes here
%   Detailed explanation goes here

b = zeros(1, params.m);
A = sparse(params.m, params.m);
 
for i = 1: params.m
    [ A , b ] = feval(func_scheme, params, A, b, i);
end

end