function [ b ] = dirichlet_boundary( b, i, alpha, beta)
% DIRICHLET_BOUNDARY 
% compute the Dirichlet boundary on the vector b function at lab(i) indice
% for m total points
% return:
%   b: the vector that contain all boundary conditions applied on the
%   system 
% parameters:
%   i: the indice with it we calculate the label
%   b: the vector that contain all boundary conditions applied on the
%   system 

if(iscolumn(b))
    s = size(b,1);
else
    s = size(b,2);
end 

if i-1 < 1
    b(i) = alpha;
end

if i+1 > s
    b(i) = beta;
end

end

