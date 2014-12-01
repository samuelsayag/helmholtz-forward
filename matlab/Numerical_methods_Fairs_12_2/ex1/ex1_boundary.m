function [ b ] = ex1_boundary( b, n, m, i, j, lab_func )
%FUNC_BOUND Summary of this function goes here
%   Detailed explanation goes here

l = feval(lab_func, i, j);

% if i-1 <= 0
%     b(l) = b(l) + 0;
% end

if j+1 >= m
    b(l) = (200 * (0.5/n) * i);
end

if i+1 >= n
    b(l) = b(l) + (200 * (0.5/m) * j);
end

% if j-1 <= 0
%     A(l, lab_func(i, j-1)) = -4;
% end

end

