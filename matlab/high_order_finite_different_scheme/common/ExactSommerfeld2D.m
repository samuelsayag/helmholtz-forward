classdef ExactSommerfeld2D
    %EXACTSOMMERFELD2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        h;
        beta;
        theta;
        scheme;
    end
    
    methods 
        function obj = ExactSommerfeld2D( h, k, theta, scheme)
        % Ord6thSommerfeld2D
        % h: the the step of the grid
        % beta: parameter of the formula (pi - k^2 for instance)
            narginchk(4, 4)
            obj = obj.check_param( h, k, theta, scheme);            
            obj.h = h;
            obj.beta = k;
            obj.theta = theta;
            obj.scheme = scheme;
        end
        
        function sx = sx( obj )
            sx =  obj.s0( @(x) x * cos(obj.theta) );   
        end
        
        function sy = sy( obj )
            sy =  obj.s0( @(x) x * sin(obj.theta) );   
        end
        
        function a0 = corner_a0(side)
            % the coefficient of the sommerfeld corner in the exact scheme
            % depend on the type of side (N= NE, E=SE, S=SW, W=NW if the
            % convention is to turn always clockwise)           
            if strcmp(side, 'north')
                a0 = - obj.scheme.a0 + obj.sx + obj.sy;
            elseif strcmp(side, 'east')
                a0 = - obj.scheme.a0 + obj.sx - obj.sy;
            elseif strcmp(side, 'south')
                a0 = - obj.scheme.a0 - obj.sx - obj.sy;
            elseif strcmp(side, 'west')
                a0 = - obj.scheme.a0 - obj.sx + obj.sy;
            else
               error('not valid side to get corner coefficient'); 
            end
                
        end
        
        function as = corner_as(obj)
            obj.h; % dummy instruction to erase warning
            as = 1;
        end        
    end
    
    methods (Access = private)     
        
        function s0 = s0( obj, f_h )
            s0 =  2 * 1i * sin(f_h(obj.beta) * obj.h);   
        end        
        
        function obj = check_param( obj, h, beta, theta, scheme )                                 
            p = inputParser;           
            
            addRequired(p, 'h', @isnumeric);   
            addRequired(p, 'beta', @isnumeric);     
            addRequired(p, 'theta', @isnumeric);
            schemes = {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D',...
                            'Ord6thHelmholtz2D', 'ExactScheme2D'};            
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'})); 
            
            parse( p, h, beta, theta, scheme);                        
        end        
    end
    
end