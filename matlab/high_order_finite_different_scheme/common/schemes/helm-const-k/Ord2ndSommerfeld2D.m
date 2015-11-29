classdef Ord2ndSommerfeld2D < SommerfeldScheme
    %ORD2NDSOMMERFELD2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        h;
        beta;
    end
    
    methods 
        function obj = Ord2ndSommerfeld2D( h, beta )
        % Ord2ndSommerfeld2D
        % h: the the step of the grid
        % beta_x: the parameter such that du_x/dx + i * beta_x = 0
        % beta_y: the parameter such that du_x/dx + i * beta_y = 0                    
            narginchk(2, 2)
            obj = obj.check_param( h, beta);            
        end
        
        function sox = sox( obj )
            sox =  obj.s0(obj.beta.x);   
        end
        
        function soy = soy( obj )
            soy =  obj.s0(obj.beta.y);   
        end 
                
        function a0 = corner_a0(side)        
            possible = {'north', 'east', 'south', 'west'};
            if strcmp(side, possible)
                bh = obj.beta * obj.h;
                a0 = -2 + 1/2 * (bh).^2 + 1i * sqrt(2) * bh; 
            end
        end
        
        function as = corner_as(obj)
            obj.h; % dummy instruction to erase warning
            as = 1;
        end         
    end
    
    methods ( Static )
        % the coefficient of the point before the central point
        function sb = sb ()
           sb = -1; 
        end 
        
        % the  coefficient of the point after the central point
        function sf = sf ()        
            sf = 1;
        end         
    end
    
    methods (Access = private)     
        
        function s0 = s0(obj, beta)
            s0 = 2 * 1i * beta * obj.h ;   
        end        
        
        function obj = check_param(obj, h, beta)                                 
            p = inputParser;           
            
            addRequired(p, 'h', @isnumeric);   
            
            function res = valide_beta(beta)
                validateattributes( beta, {'struct'}, {'nonempty'});
                resx = isfield(beta, 'x');
                if resx
                    validateattributes( beta.x, {'numeric'}, {'nonempty'});
                end
                resy = isfield(beta, 'y');
                if resy
                    validateattributes( beta.y, {'numeric'}, {'nonempty'});                
                end
                res = resx || resy;
            end
            addRequired(p, 'beta', @(x)valide_beta(x));     
            
            parse( p, h, beta );            
            
            obj.h = h;
            obj.beta = beta;
        end        
    end
    
end

