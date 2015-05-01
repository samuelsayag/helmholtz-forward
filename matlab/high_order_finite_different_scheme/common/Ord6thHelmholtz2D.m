classdef Ord6thHelmholtz2D
    %ORD6THHELMHOLTZ2D this class contain the 4th order scheme for the
    %Helmholtz equation.
    %
    %Reference:
    %2012 Erlangga, Turkel - ITERATIVE SCHEMES FOR HIGH ORDER COMPACT 
    %DISCRETIZATIONS
    
    properties (SetAccess = public)
        k = 0;
        h = 0;
        delta = 0;
    end
    
    methods (Access = public)
        
        function obj = Ord6thHelmholtz2D(k, h, d)
        % Ord2ndHelmholtz
        % k: the k of an Helmholtz equation 
        % h: the step length (basic division of the grid)
        % delta: a parameter of the scheme
            narginchk(2, 3)
            if nargin == 2
                d = 0;
            end
            obj = obj.check_param(k, h, d);
        end
        
        function a0 = a0(obj)
        % return A0 coefficient
            kh = obj.k * obj.h;            
            a0 = -10/3 + kh^2 * 67/90 + kh^4 * (obj.delta - 3)/180; 
        end
        
        function as = as(obj)
        % return the As coefficient
            kh = obj.k * obj.h;    
            as = 2/3 + kh^2 * 2/45 + kh^4 * (3 - 2 * obj.delta)/720; 
        end
        
        function ac = ac(obj)
        % return the Ac coefficient
            kh = obj.k * obj.h;
            ac = 1/6 + kh^2 * 7/360 + kh^4 * obj.delta/144; 
        end

        function obj = check_param(obj, k, h, d)            
            p = inputParser;
            addRequired(p, 'k', @isnumeric);
            addRequired(p, 'h', @isnumeric);
            addRequired(p, 'd', @isnumeric);
            parse(p, k, h, d);            
            obj.k = k;
            obj.h = h;
            obj.delta = d;                                          
        end
    end    
end

