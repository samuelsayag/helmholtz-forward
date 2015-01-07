function [ A , b ] = dirichlet_sommerfeld_boundary( A, b, i, k, h, alpha,...
    left_dirichlet)
% DIRICHLET_SOMMERFELD_BOUNDARY 
% compute the Dirichlet and the Sommerfeld boundary on the vector b 
% function at lab(i) indice for m total points
% return:
%   b: the vector that contain all boundary conditions applied on the
%   system 
% parameters:
%   i: the indice with it we calculate the label
%   b: the vector that contain all boundary conditions applied on the
%   system 

m = size(A,2);

%-----------------------------------------
% INTERIOR POINTS
%-----------------------------------------
if i > 1 && i < m 
    A(i, i) = 2 - (k * h).^2;

    if i-1 > 0
        A(i, i-1) = -1;
    end

    if i+1 <= m
        A(i, i+1) = -1;
    end
end


%-----------------------------------------
% BOUNDARY POINTS
% DIRICHLET AND SOMMERFELD
%-----------------------------------------

if(iscolumn(b))
    s = size(b,1);
else
    s = size(b,2);
end 

% left dirichlet and right sommerfeld
if left_dirichlet        
    if i == 1
        A(i, i) = 2 - (k * h).^2;
        A(i, i+1) = -1;
        b(i) = alpha;
    end
    if i == s
        A(i, i) = 2 - (k * h).^2 - 1i * 2 * sin(k*h);
        A(i, i-1) = -2;        
    end
else % right dirichlet and left sommerfeld       
    if i == s
        A(i, i) = 2 - (k * h).^2;
        A(i, i-1) = -1;
        b(i) = alpha;
    end
    if i == 1
        A(i, i) = 2 - (k * h).^2 + 1i * 2 * sin(k*h);
        A(i, i+1) = -2;        
    end
end 

end