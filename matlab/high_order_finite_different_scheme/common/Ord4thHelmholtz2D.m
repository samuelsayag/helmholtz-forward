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
        
        function obj = Ord4thHelmholtz2D(k, h, g)
        % Ord2ndHelmholtz
        % k: the k of an Helmholtz equation 
        % h: the step length (basic division of the grid)
        % gamma: a parameter of the scheme
            narginchk(2, 3)
            if nargin == 2
                g = 0;
            end
            obj = obj.check_param(k, h, g);
        end
        
        function a0 = a0(obj)
            % return A0 coefficient
            a0 = -10/3 + (obj.k * obj.h)^2 * (2/3 + obj.gamma/36); 
        end
        
        function as = as(obj)
            % return the As coefficient
            as = 2/3 + (obj.k * obj.h)^2 * (1/12 - obj.gamma/72); 
        end
        
        function ac = ac(obj)
           % return the Ac coefficient
           ac = 1/6 + (obj.k * obj.h)^2 * (obj.gamma/144); 
        end
        
        function bs = bs(obj)
            % return the coefficient of the dirichlet point on SIDE
            obj.h; % dummy instruction 
            bs = 1; 
        end
                
        function bc = bc(obj)
            % return the coefficient of the dirichlet point on CORNER
            obj.h; % dummy instruction
            bc = 0; 
        end
        
        function obj = check_param(obj, k, h, g)            
            p = inputParser;
            addRequired(p, 'k', @isnumeric);
            addRequired(p, 'h', @isnumeric);
            addRequired(p, 'g', @isnumeric);
            parse(p, k, h, g);            
            obj.k = k;
            obj.h = h;
            obj.gamma = g;                                          
        end
    end    
end

