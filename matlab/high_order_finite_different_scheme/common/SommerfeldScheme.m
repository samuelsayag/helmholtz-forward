classdef (Abstract) SommerfeldScheme
    %CENTRALSCHEME abstract class that specifies the methods of a central
    %scheme for the finite difference framework.
    
    methods (Abstract)
        % the horizontal side point coefficient
        sx = sx( obj );
        % the vertical side point coefficient
        sy = sy( obj );
        % the central corner point coefficient
        a0 = corner_a0(obj, side);
        % the side corner point coefficient
        as = corner_as(obj);
    end
    
end

