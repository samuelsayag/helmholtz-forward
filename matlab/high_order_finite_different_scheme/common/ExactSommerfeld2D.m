classdef ExactSommerfeld2D
    %EXACTSOMMERFELD2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        h;
        beta;
        theta;
    end
    
    methods 
        function obj = ExactSommerfeld2D( h, k, theta )
        % Ord6thSommerfeld2D
        % h: the the step of the grid
        % beta: parameter of the formula (pi - k^2 for instance)
            narginchk(3, 3)
            obj = obj.check_param( h, k, theta);            
        end
        
        function sx = sx( obj )
            sx =  obj.s0( @(x) x * cos(obj.theta) );   
        end
        
        function sy = sy( obj )
            sy =  obj.s0( @(x) x * sin(obj.theta) );   
        end        
        
    end
    
    methods (Access = private)     
        
        function s0 = s0( obj, f_h )
            s0 =  2 * 1i * sin(f_h(obj.beta) * obj.h);   
        end        
        
        function obj = check_param( obj, h, beta, theta )                                 
            p = inputParser;           
            
            addRequired(p, 'h', @isnumeric);   
            addRequired(p, 'beta', @isnumeric);     
            addRequired(p, 'theta', @isnumeric);
            
            parse( p, h, beta, theta);            
            
            obj.h = h;
            obj.beta = beta;
            obj.theta = theta;
        end        
    end
    
end

