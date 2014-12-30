function [ b ] = dirichlet_boundary( b, n, m, i, j, lab_func )
%FUNC_BOUND Summary of this function goes here
%   Detailed explanation goes here

l = feval(lab_func, i, j);
b(l) = 0;


 if i-1 <= 0
     b(l) = b(l) + 0;
 end

if j+1 >= m
    b(l) = b(l) + (200 * (0.5/n) * i);
end

if i+1 >= n
    b(l) = b(l) + (200 * (0.5/m) * j);
end

if j-1 <= 0
    b(l) = b(l) + 0;
end

end

