function [ A, b ] = build_two_dimensional_problem2( param, func_scheme )
%BUILD_SQUARE_GRID Summary of this function goes here
%   Detailed explanation goes here

% m: nbr of column of the domain, n: nbr of line of the domain
L = param.n * param.m;
% matrix A and b are declared sparse for memory space saving reasons.
b = sparse(L, 1);
A = sparse(L, L);

for i = 1: param.m
    for j = 1:param.n
        % The function pointer func_scheme is called withe the paramters: 
        % param(structure that contains the simulation parameter), A, b, i,
        % and j. A and b are progressively "enrich" with a new line for 
        % each couple of index (i,j).
        [ A, b ] = feval(func_scheme, param, A, b, i , j);
    end
end
   
end
