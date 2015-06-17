classdef (Abstract) CentralScheme
    %CENTRALSCHEME abstract class that specifies the methods of a central
    %scheme for the finite difference framework.
    
    methods (Abstract)
        % the central point coefficient
        a0 = a0(obj);
        % the side point coefficient
        as = as(obj);
        % the corner point coefficient
        ac = ac(obj);
    end
    
end

