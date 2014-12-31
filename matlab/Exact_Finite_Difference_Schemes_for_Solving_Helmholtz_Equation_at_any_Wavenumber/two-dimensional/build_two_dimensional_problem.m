function [ A, b ] = build_two_dimensional_problem( n, m, func_inter, func_bound )
%BUILD_SQUARE_GRID Summary of this function goes here
%   Detailed explanation goes here


L = n * m;
b = zeros(1, L);
A = sparse(L, L);
 
lab = @(i,j) i + (m-j) * n; 

for i = 1:n    
    for j = 1:m
        A = feval(func_inter, A, n, m, i , j, lab);
        b = feval(func_bound, b, n, m, i, j, lab);
    end
end
   
end

