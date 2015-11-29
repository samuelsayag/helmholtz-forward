classdef Order6thVarK_F < F
    %ORDER6THVARK_F See the documentation of the class it inherits for
    %basic knowledge of this class.
    %This particular one is an implementation of the Abstract class F and
    %provide a scheme for the function F (when we build AU=F)in relation to
    %the scheme Ord6thVarKHelmholtz2D.
    
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
        F;
        Fx; 
        Fy; 
        Fxx; 
        Fyy;
        Fxxxx;
        Fyyyy;
        Fxxyy;
    end
    
    methods (Access = public)
        
        function obj = Order6thVarK_F( h, F, Fx, Fy, Fxx, Fyy, ...
                Fxxxx, Fyyyy, Fxxyy, K2x, K2y)        
            narginchk(9, 9);
            [obj.h, obj.F, obj.Fx, obj.Fy, ...
                obj.Fxx, obj.Fyy, obj.Fxxxx, obj.Fyyyy, obj.Fxxyy ] = ...
                obj.check_param2( h, F, Fx, Fy, Fxx, Fyy, ...
                Fxxxx, Fyyyy, Fxxyy, K2x, K2y );
        end
        
        
        % if the scheme need to be aware of the position of the indexes i
        % and j in the grid. So as to compute an (i,j) aware value
        % function (ex: k(x,y), kx(x,y)...).
        function obj = set_pos(obj, i, j)
            obj.i = i;
            obj.j = j;
        end          
        
        % provide the value of the function given a position in the grid
        function val = f(obj)
            x = obj.i;
            y = obj.j;
            val = obj.h.^2 * (...
                ...
                (1 - (obj.k2 * obj.h.^2)/20) * obj.F(x, y) + ...
                ...
                ((obj.h.^2)./12) * (obj.Fxx(x,y) + obj.Fyy(x,y)) + ...
                ...
                ((obj.h.^4)./360) * (obj.Fxxxx(x,y) + obj.Fyyyy(x,y)) + ...
                ...
                ((obj.h.^4)./90) * obj.Fxxyy(x,y) + ...
                ...
                ((obj.h.^4)./60) * ( obj.K2x(x,y) * obj.Fx(x,y) + ...
                                    obj.K2x(x,y) * obj.Fx(x,y) ));
        end
    end
    
    methods (Static, Access = private)
        function [h, k] = check_param1(h, k)            
            p = inputParser;
            addRequired(p, 'h', @isnumeric);
            addRequired(p, 'k', @ischar)                        
            parse(p, h, k);
        end        
        
        function [h, F, Fx, Fy, Fxx, Fyy, Fxxxx, Fyyyy, Fxxyy, K2x, K2y] ...
                = check_param2( h, F, Fx, Fy, Fxx, Fyy, ...
                Fxxxx, Fyyyy, Fxxyy, K2x, K2y )            
            
            p = inputParser;
            
            addRequired(p, 'h', @isnumeric);
            function r = valid_handle(x)
                validateattributes( x, {'function_handle'}, {'nonempty'});
                r = nargin(x) == 2;
            end
            addRequired(p, 'Fx', @(x)valid_handle(x))            
            addRequired(p, 'Fxy', @(x)valid_handle(x))            
            addRequired(p, 'Fxx', @(x)valid_handle(x))            
            addRequired(p, 'Fyy', @(x)valid_handle(x))            
            addRequired(p, 'Fxxxx', @(x)valid_handle(x))            
            addRequired(p, 'Fyyyy', @(x)valid_handle(x))            
            addRequired(p, 'Fxxyy', @(x)valid_handle(x))            
            addRequired(p, 'K2x', @(x)valid_handle(x))            
            addRequired(p, 'K2y', @(x)valid_handle(x))            
                        
            parse(p, h, F, Fx, Fy, Fxx, Fyy, ...
                Fxxxx, Fyyyy, Fxxyy, K2x, K2y);
        end        
    end 
    
    
end

