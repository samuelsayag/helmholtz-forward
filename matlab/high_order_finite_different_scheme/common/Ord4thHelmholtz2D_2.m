classdef Ord4thHelmholtz2D_2  < CentralScheme
    %ORD6THHELMHOLTZ2D this class contain the 4th order scheme for the
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
        
        function obj = Ord4thHelmholtz2D_2(k, h)
        % Ord2ndHelmholtz
        % k: the k of an Helmholtz equation 
        % h: the step length (basic division of the grid)
        % delta: a parameter of the scheme
            narginchk(2, 2)
            obj = obj.check_param(k, h);
        end
        
        function a0 = a0(obj)
        % return A0 coefficient
            kh = obj.k * obj.h;            
%             a0 = -10/3 + kh^2 * 41/45 - kh^4 * 1/20; 
            a0 = -10/3 + kh^2 * 41/45;
        end
        
        function as = as(obj)
        % return the As coefficient
            kh = obj.k * obj.h;    
            as = 2/3 + kh^2 * 1/90; 
        end
        
        function ac = ac(obj)
        % return the Ac coefficient
            kh = obj.k * obj.h;
            ac = 1/6 + kh^2 * 1/90; 
        end
    end
    
    methods (Access = private)
        function obj = check_param(obj, k, h)            
            p = inputParser;
            addRequired(p, 'k', @isnumeric);
            addRequired(p, 'h', @isnumeric);
            parse(p, k, h);            
            obj.k = k;
            obj.h = h;
        end    
    end
    
end

