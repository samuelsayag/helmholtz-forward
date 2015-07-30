classdef Ord6thSommerfeld2D < SommerfeldScheme
    %ORD6THHELMHOLTZSOMMERFELD2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        h;
        beta;
    end
    
    methods 
        function obj = Ord6thSommerfeld2D( h, beta )
        % Ord6thSommerfeld2D
        % h: the the step of the grid
        % beta: parameter of the formula (pi - k^2 for instance)
            narginchk(2, 2)
            obj = obj.check_param( h, beta );            
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
                a0 = 0; 
            end
        end
        
        function as = corner_as(obj)
            obj.h; % dummy instruction to erase warning
            as = 1;
        end         
        
    end
    
    methods (Static, Access = public )
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
            s0 = 2 * 1i * beta * obj.h * ( 1 ...
                - (beta * obj.h).^2/6 ...
                + (beta * obj.h).^4/120 );  
        end        
        
        function obj = check_param(obj, h, beta )                                 
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