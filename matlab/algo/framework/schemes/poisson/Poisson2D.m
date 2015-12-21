classdef Poisson2D  < CentralScheme
    %POISSON2D Recoding of the Poisson equation with the new framework     
    methods (Static)
        % the central point coefficient
        function a0 = a0()
           a0 = 4; 
        end
        % the side point coefficient
        function as = as()
           as = -1; 
        end
        % the corner point coefficient
        function ac = ac()
            ac = 0;
        end 
    end    
end

