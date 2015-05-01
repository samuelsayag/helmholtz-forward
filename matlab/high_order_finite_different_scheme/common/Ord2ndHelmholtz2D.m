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
        
        function obj = Ord2ndHelmholtz2D(k, h)
        % Ord2ndHelmholtz
        % k = the k of an Helmholtz equation 
        % h = the step length (basic division of the grid)            
            obj = obj.check_param(k, h);
        end        
        
        function a0 = a0(obj)
        % return A0 coefficient
            a0 = - 4 + (obj.k * obj.h)^2; 
        end
        
        function as = as(obj)
        % return the As coefficient
            obj.h; % dummy instruction    
            as = 1; 
        end
        
        function ac = ac(obj)
        % return the Ac coefficient
            obj.h; % dummy instruction
            ac = 0; 
        end
        
        function obj = check_param(obj, k, h)
            narginchk(3, 3)
            p = inputParser;
            addRequired(p, 'k', @isnumeric);
            addRequired(p, 'h', @isnumeric);        
            parse(p, k, h);
            obj.k = k;
            obj.h = h;            
        end
    end
    
end

