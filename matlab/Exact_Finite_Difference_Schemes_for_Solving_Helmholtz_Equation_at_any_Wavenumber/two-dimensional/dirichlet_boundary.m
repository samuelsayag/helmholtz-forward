function [ b ] = dirichlet_boundary( b, n, m, i, j, lab_func, S_dirchlet,... 
W_dirichlet, N_dirichlet, E_dirichlet)
% DIRICHLET_BOUNDARY 
% compute the Dirichlet boundary on the vector b function at lab(i) indice
% for m total points
% return:
%   b: the vector that contain all boundary conditions applied on the
%   system 
% parameters:
%   b: the vector that contain all boundary conditions applied on the
%   system
%   n, m: number of line of the matrix, number of column of the system
%   matrix
%   i, j: the indices of respectively line and column in the system matrix
%   lab_func: create the label in the b vector from the point
%   X_dirichlet: X may be (Sout,West,North,East) and is the dirichlet
%   imposed on the side of the system.

l = feval(lab_func, i, j);
b(l) = 0;

% South boundary
if j-1 < 1
    b(l) = b(l) + S_dirchlet;
end
% West boundary
if i-1 < 1
    b(l) = b(l) + W_dirichlet;
end
% North boundary
if j+1 > m
    b(l) = b(l) + N_dirichlet;
end
% East boundary
if i+1 > n
    b(l) = b(l) + E_dirichlet;
end

end

