classdef F
    %F This is the general interface for a scheme of a function F such that
    %the problem to solve is: (Opertor)U(x,y) = F. For instance in the
    %Helmholtz equation: D²u + k²u = F(x,y) (D² is the Laplacian operator).
    %With finite differenc it gives birth to a linear system of
    %equation: AU = F (U, a vector, is the search solution, A, a n-diagonal
    %matrix generally is the result of the scheme used on all the grid of
    %the model and F is vector of the value function at each point.
    
    methods (Abstract)
        % if the scheme need to be aware of the position of the indexes i
        % and j in the grid. So as to compute an (i,j) aware value
        % function (ex: k(x,y), kx(x,y)...).
        set_pos(i,j);        
        % provide the value of the function given a position in the grid
        f = val(obj)
    end
    
end

