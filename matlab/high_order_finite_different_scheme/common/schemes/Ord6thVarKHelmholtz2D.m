classdef Ord6thVarKHelmholtz2D  < NinePtStencil
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
        str_k;
        k;
        Kx; 
        Ky; 
        Kxx; 
        Kyy;
    end
    
    methods (Access = public)
        
        function obj = Ord6thVarKHelmholtz2D( h, k, Kx, Ky, Kxx, Kyy)        
            narginchk(6, 6);
            [ obj.k, obj.Kx, obj.Ky, obj.Kxx, obj.Kyy ] = ...
                check_param2( h, k, Kx, Ky, Kxx, Kyy );
        end
        
        % if the scheme need to be aware of the position of the indexes i
        % and j in the grid. So as to compute an (i,j) aware value
        % function (ex: k(x,y), kx(x,y)...).
        function set_pos(obj, i, j)
            obj.i = i;
            obj.j = j;
        end             
        
        % the central point coefficient
        function c = c(obj)
        end
        % the north point coefficient
        function n = n(obj)
        end
        % the north east point coefficient
        function ne = ne(obj)
        end
        % the east point point coefficient
        function e = e(obj)
        end
        % the south east point coefficient
        function se = se(obj)
        end
        % the south point coefficient
        function s = s(obj)
        end
        % the south west point coefficient
        function sw = sw(obj)
        end
        % the west point coefficient
        function w = w(obj)
        end
        % the north west point coefficient
        function nw = nw(obj)
        end
        
    end
    
    methods (Access = private)
        
        function a0 = a0(obj)
            kh = obj.k * obj.h;
            D2 = (obj.Kxx(obj.i,obj.j) + obj.Kyy(obj.i,obj.j));
            a0 = -10/3 ...
                + kh^2 * 41/45 ...
                - kh^4 * 1/20 ...
                + D2 .* obj.h.^4/20; 
        end
        
        function as = as(obj)
            as = 2/3 + (obj.k * obj.h).^2 * 1/90; 
        end
        
        function ac = ac(obj)
            ac = 1/6 + (obj.k * obj.h).^2 * 1/90; 
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
            addRequired(p, 'Kx', @(x)valid_handle(x))            
            addRequired(p, 'Ky', @(x)valid_handle(x))            
            addRequired(p, 'Kxx', @(x)valid_handle(x))            
            addRequired(p, 'Kyy', @(x)valid_handle(x))            
                        
            parse(p, h, k, Kx, Ky, Kxx, Kyy);
        end        
    end 
    
    methods (Static, Access = public)
        function [ k, Kx, Ky, Kxx, Kyy ] = ...
                build_derivative(k_str)
            
            syms x y; % declare 2 symbolic variables
            toHandle = @(x) str2func(strcat('@(x,y)',char(x)));
            
            ks = sym(k_str);
            k = toHandle(ks);
            
            kxs = diff( ks, x );
            Kx = toHandle(kxs);
            
            kys = diff( ks, y );
            Ky = toHandle(kys);
            
            kxxs = diff(kxs);
            Kxx = toHandle(kxxs);
            
            kyys = diff(kys);
            Kyy = toHandle(kyys);            
        end        
    end
    
end

