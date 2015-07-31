classdef Ord2ndVarKHelmholtz2D  < NinePtStencil
    %ORD6THHELMHOLTZ2D this class contain the 4th order scheme for the
    %Helmholtz equation.
    %
    %Reference:
    %2012 Erlangga, Turkel - ITERATIVE SCHEMES FOR HIGH ORDER COMPACT 
    %DISCRETIZATIONS
    
    properties (SetAccess = private)
        % the position of the stencil in the grid. Is of use to calculate
        % postion aware function such as k,kx...
        i; % column index (the 'x' axis index)
        j; % line index (the 'y' axix index)
        
        % Ord2ndHelmholtz
        % (1) K, (2) Kx, (3) Ky, (4) Kxx, (5) Kyy : 
        % respectively (1) a function of (potentially) variables x and y 
        % (2D), (2)(3) its first derivative with respect to x and y and (3)
        % and (4) its second derivative with respect to x and y.
        % h: the step length (basic division of the grid)
        h;
        k2;
    end
    
    methods (Access = public)
        
        function obj = Ord6thVarKHelmholtz2D( h, k2)        
            narginchk(6, 6);
            [obj.h, obj.k2] = obj.check_param2( h, k2);
        end
        
        % if the scheme need to be aware of the position of the indexes i
        % and j in the grid. So as to compute an (i,j) aware value
        % function (ex: k(x,y), kx(x,y)...).
        function obj = set_pos(obj, i, j)
            obj.i = i;
            obj.j = j;
        end             
        
        % the central point coefficient
        function c = c(obj)
            c = obj.a0;
        end
        % the north point coefficient
        function n = n(obj)
            n = obj.as + obj.bs('north');
        end
        % the north east point coefficient
        function ne = ne(obj)
            ne = obj.ac + obj.bc(+1, +1);
        end
        % the east point point coefficient
        function e = e(obj)
            e = obj.as + obj.bs('east');
        end
        % the south east point coefficient
        function se = se(obj)
            se = obj.ac + obj.bc(+1, -1);
        end
        % the south point coefficient
        function s = s(obj)
            s = obj.as + obj.bs('south');
        end
        % the south west point coefficient
        function sw = sw(obj)
            sw = obj.ac + obj.bc(-1, -1);
        end
        % the west point coefficient
        function w = w(obj)
            w = obj.as + obj.bs('west');
        end
        % the north west point coefficient
        function nw = nw(obj)
            nw = obj.ac + obj.bc(-1, +1);
        end
        
    end
    
    methods (Access = private)
        
        function a0 = a0(obj)
            kh2 = obj.k2(obj.i,obj.j) * obj.h.^2;
            D2 = (obj.K2xx(obj.i,obj.j) + obj.K2yy(obj.i,obj.j));
            a0 = -10/3 ...
                + kh2 * 41/45 ...
                - kh2^2 * 1/20 ...
                + D2 .* obj.h.^4/20; 
        end
        
        function as = as(obj)
            kh2 = obj.k2(obj.i,obj.j) * obj.h.^2;
            as = 2/3 + kh2 * 1/90; 
        end
        
        function ac = ac(obj)
            kh2 = obj.k2(obj.i,obj.j) * obj.h.^2;
            ac = 1/6 + kh2 * 1/90; 
        end

        function bc  = bc(obj, s1, s2)
            dk = s1 * obj.K2x(obj.i,obj.j) + s2 * obj.K2y(obj.i,obj.j);
            bc = dk * obj.h.^3/120;
        end
        
        function bs  = bs(obj, axis)
            if strcmp(axis,'north')
                dk = obj.K2y(obj.i,obj.j);
            elseif strcmp(axis,'south')
                dk = - obj.K2y(obj.i,obj.j);
            elseif strcmp(axis,'east')
                dk = obj.K2x(obj.i,obj.j);
            elseif strcmp(axis,'west')
                dk = -obj.K2x(obj.i,obj.j);
            end
            
            ck = 2/3 + (obj.k2(obj.i, obj.j) .* obj.h^2)/6;            
            bs = obj.h.^3/20 .* ck * dk;
        end        

    end
    
    methods (Static, Access = private)
        function [h, k] = check_param1(h, k)            
            p = inputParser;
            addRequired(p, 'h', @isnumeric);
            addRequired(p, 'k', @ischar)                        
            parse(p, h, k);
        end        
        
        function [h, k, Kx, Ky, Kxx, Kyy] ...
                = check_param2( h, k, Kx, Ky, Kxx, Kyy )            
            p = inputParser;
            
            addRequired(p, 'h', @isnumeric);
            function r = valid_handle(x)
                validateattributes( x, {'function_handle'}, {'nonempty'});
                r = nargin(x) == 2;
            end
            addRequired(p, 'k', @(x)valid_handle(x))                    
                        
            parse(p, h, k);
        end        
    end 
    
    methods (Static, Access = public)
        function [k2] = build_derivative(k, param)
            % build the different derivative of k² needed for this scheme
            % to comute the coefficient.
            syms x y; % declare 2 symbolic variables
            tof = @( z ) matlabFunction( z , 'Vars', [x,y]);            
            k2s = sym(k).^2;
            k2 = tof(k2s);                  
            
            k2 = discreteWrapper(param.a, param.c, param.h, k2);
            
        end
    end
    
end

