function [ imp_x ] = improved_x( A, b, x, omega )
% Gives an improved vector imp_x from a x vector and a relaxation term.
% USAGE: [ imp_x ] = improved_x( x, omega )
% INPUT:
% x = starting solution vector
% omega = computed relaxation factor
% OUTPUT:
% A = A square matrix (n * n)
% b = a vector (1 * n) or (n * 1)
% x = solution vector
% omega = relaxation factor
flip = 0;
if ~iscolumn(x)
    x = x';
    flip = 1;
end

if ~iscolumn(b)
    b = b';    
end

d = diag(A);
s = (A * x) - (d .* x);
imp_x = (1 ./ d) .* (b - s);  
imp_x = omega .* imp_x + (1 - omega) * x;  

if flip
    imp_x = imp_x';
else
    imp_x = imp_x;
end

end

