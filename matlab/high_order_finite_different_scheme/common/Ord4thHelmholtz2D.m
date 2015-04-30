classdef Ord4thHelmholtz2D
    %ORD4THHELMHOLTZ2D this class contain the 4th order scheme for the
    %Helmholtz equation.
    %
    %Reference:
    %2012 Erlangga, Turkel - ITERATIVE SCHEMES FOR HIGH ORDER COMPACT 
    %DISCRETIZATIONS
    
    properties (SetAccess = public)
        k = 0;
        h = 0;
        gamma = 0;
    end
    
    methods (Access = public)
        
        function obj = Ord2ndHelmholtz(k, h, g)
        % Ord2ndHelmholtz
        % k: the k of an Helmholtz equation 
        % h: the step length (basic division of the grid)
        % gamma: a parameter of the scheme 
            if nargin == 0
                obj.k = 0;
                obj.h = 0;
                obj.gamma = 0;
            elseif nargin == 1
                obj.k = k;
                obj.h = 0;
                obj.gamma = 0;
            elseif nargin == 2
                obj.k = k;
                obj.h = h;
                obj.gamma = 0;
            else
                obj.k = k;
                obj.h = h;
                obj.gamma = g;                
            end
        end
        
        function a0 = a0(obj)
        % return A0 coefficient
            a0 = -10/3 + (obj.k * obj.h)^2 * (2/3 + obj.gamma/36); 
        end
        
        function as = as(obj)
        % return the As coefficient
            as = 2/3+ (obj.k * obj.h)^2 * (1/12 - obj.gamma/72); 
        end
        
        function ac = ac(obj)
        % return the Ac coefficient
           ac = 1/6 + (obj.k * obj.h)^2 * (obj.gamma/144); 
        end
        
    end    
end

