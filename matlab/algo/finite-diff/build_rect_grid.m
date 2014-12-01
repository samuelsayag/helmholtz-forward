function [ A, b ] = build_rect_grid( n, m, func_inter, func_bound )
%BUILD_SQUARE_GRID Summary of this function goes here
%   Detailed explanation goes here

n = n+1;
m = m+1;

L = (n-1) * (m-1);
b = zeros(1, L);
A = sparse(L, L);
 
lab = @(i,j) i+(m-1-j)*(n-1); 

for i = 1:n-1
    for j = 1:m-1
        A = feval(func_inter, A, n, m, i , j, lab);
        b = feval(func_bound, b, n, m, i, j, lab);
    end
end
   
end