classdef Ord2ndHelmholtz2D
    %ORD2NDHELMHOLTZ2D this class contain the 2nd order scheme for the
    %Helmholtz equation.
    %
    %Reference:
    %2012 Erlangga, Turkel - ITERATIVE SCHEMES FOR HIGH ORDER COMPACT 
    %DISCRETIZATIONS    
    
    properties (SetAccess = public)
        k = 0;
        h = 0;
    end
    
    methods (Access = public)
        
        function obj = Ord2ndHelmholtz(k, h)
        % Ord2ndHelmholtz
        % k = the k of an Helmholtz equation 
        % h = the step length (basic division of the grid)
            if nargin == 0
                obj.k = 0;
                obj.h = 0;
            elseif nargin == 1
                obj.k = k;
                obj.h = 0;
            else
                obj.k = k;
                obj.h = h;
            end
        end        
        
        function a0 = a0(obj)
        % return A0 coefficient
            a0 = - 4 + (obj.k * obj.h)^2; 
        end
        
        function as = as(obj)
        % return the As coefficient
            as = (obj.k * obj.h)^2 - 4; 
        end
        
        function ac = ac(obj)
        % return the Ac coefficient
            obj.h; % dummy instruction
            ac = 0; 
        end
        
    end
    
end

