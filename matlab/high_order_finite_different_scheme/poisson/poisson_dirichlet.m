function [ value ] = poisson_dirichlet( x, y )
%POISSON_DIRICHLET definition of Dirichlet Boundary Condition
    value = 0;
    if x < 0.125 
        value = 0;
    end
    
    if x > 0.375
        value = value + 200 * y;
    end 
    
    if y < 0.125
        value = value + 0;
    end 

    if y > 0.375 
        value = value + 200 * x;
    end
end

