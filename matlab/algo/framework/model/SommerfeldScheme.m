classdef (Abstract) SommerfeldScheme
    %CENTRALSCHEME abstract class that specifies the methods of a central
    %scheme for the finite difference framework.
    
    methods (Abstract)
        % the horizontal side point coefficient
        sox = sox( obj );
        % the vertical side point coefficient
        soy = soy( obj );
        % the point before the central point
        sb = sb ( obj );
        % the point after the central point
        sf = sf ( obj );
        % the central corner point coefficient
        a0 = corner_a0(obj, side);
        % the side corner point coefficient
        as = corner_as(obj);
    end
    
end

