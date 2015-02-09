function [ A, b ] = build_two_dimensional_problem2( params, func_scheme )
%BUILD_SQUARE_GRID Summary of this function goes here
%   Detailed explanation goes here

L = params.n * params.m;
b = sparse(1, L);
A = sparse(L, L);

for i = 1: params.m    
    for j = 1:params.n
        [A,b] = feval(func_scheme, params, A, b, i , j);
    end
end
   
end
