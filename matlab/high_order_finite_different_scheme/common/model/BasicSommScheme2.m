classdef BasicSommScheme2
    %BASICSOMMSCHEME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        scheme;
        sommerfeld;
    end
    
    methods ( Access = public )
        
        function obj = BasicSommScheme( scheme, sommerfeld )
        % scheme : basic scheme BasicScheme (loaded with a specific central
        % scheme.
        % sommerfeld: the sommerfeld scheme for the boundary condition
            narginchk(2, 2)
            obj = obj.check_param( scheme, sommerfeld );            
        end        
        
        function v_A = n_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % east point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south east point
            v_A(4) = obj.scheme.as .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south point
            v_A(5) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south west point
            v_A(6) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % west point            
        end
        
        function v_A = e_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % central point
            v_A(2) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % north point
            v_A(3) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % south point
            v_A(4) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south west point
            v_A(5) = obj.scheme.as .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % west point            
            v_A(6) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % north west point
        end
        
        function v_A = s_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.as .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north east point
            v_A(4) = obj.scheme.as -...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % east point
            v_A(5) = obj.scheme.as -...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % west point
            v_A(6) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north west point
        end        
        
        function v_A = w_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % north point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north east point
            v_A(4) = obj.scheme.as .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % east point
            v_A(5) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % south east point
            v_A(6) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % south point
        end        
        
        function v_A = n_half_ne_pt( obj )             
            v_A = zeros(4,1);                                    
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.as .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south west point
            v_A(4) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % west point
        end
        
        function v_A = e_half_ne_pt( obj )   
            v_A = zeros(4,1);                                    
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % south point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south west point
            v_A(4) = obj.scheme.as .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % west point
        end
        
        function v_A = e_half_se_pt( obj )
            v_A = [];             
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % north point
            v_A(3) = obj.scheme.as .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % west point
            v_A(4) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % north west point            
        end
        
        function v_A = s_half_se_pt( obj )            
            v_A = [];             
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % central point            
            v_A(2) = obj.scheme.as .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north point
            v_A(3) = obj.scheme.as -...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % west point
            v_A(4) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north west point
        end
        
        function v_A = s_half_sw_pt( obj )
            v_A = [];
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % central point            
            v_A(2) = obj.scheme.as .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north east point
            v_A(4) = obj.scheme.as -...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % east point            
        end        
        
        function v_A = w_half_sw_pt( obj )            
            v_A = [];
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % central point            
            v_A(2) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % north point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % north east point
            v_A(4) = obj.scheme.as .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % east point                        
        end        
        
        function v_A = w_half_nw_pt( obj )            
            v_A = [];
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % central point            
            v_A(2) = obj.scheme.as .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % east point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sf ./ obj.sommerfeld.sb); % south east point
            v_A(4) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % south point            
        end         
        
        function v_A = n_half_nw_pt( obj )            
            v_A = [];
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.as - ...
                obj.scheme.ac .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % east point
            v_A(3) = obj.scheme.ac .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south east point
            v_A(4) = obj.scheme.as .* (1 - obj.sommerfeld.sb ./ obj.sommerfeld.sf); % south point            
        end
        
        function v_A = ne_pt(obj)
            v_A = zeros(4,1);
            v_A(1) = obj.sommerfeld.corner_a0('north'); % central point            
            v_A(2) = obj.sommerfeld.corner_as(); % south point
            v_A(3) = 0; % south west point
            v_A(4) = obj.sommerfeld.corner_as(); % west point    
        end

        function v_A = se_pt(obj)
            v_A = zeros(4,1);
            v_A(1) = obj.sommerfeld.corner_a0('east'); % central point            
            v_A(2) = obj.sommerfeld.corner_as(); % north point
            v_A(3) = obj.sommerfeld.corner_as(); % west point
            v_A(4) = 0; % north west point    
        end

        function v_A = sw_pt(obj)
            v_A = zeros(4,1);                                    
            v_A(1) = obj.sommerfeld.corner_a0('south'); % central point            
            v_A(2) = obj.sommerfeld.corner_as(); % north point
            v_A(3) = 0; % north east point
            v_A(4) = obj.sommerfeld.corner_as(); % east point    
        end

        function v_A = nw_pt(obj)
            v_A = zeros(4,1);
            v_A(1) = obj.sommerfeld.corner_a0('west'); % central point            
            v_A(2) = obj.sommerfeld.corner_as(); % east point
            v_A(3) = 0; % south east point
            v_A(4) = obj.sommerfeld.corner_as(); % south point    
        end            
    end

    methods ( Access = private )
        
        function obj = check_param(obj, scheme, sommerfeld)                                 
            p = inputParser;

            schemes = {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D',...
                'Ord6thHelmholtz2D', 'ExactScheme2D'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));            
            
            sommerfelds = {'Ord6thSommerfeld2D', 'Ord2ndSommerfeld2D', ...
                'ExactSommerfeld2D'};
            addRequired(p, 'sommerfeld', ...
                @(x)validateattributes( x, sommerfelds, {'nonempty'}));
            
            parse(p, scheme, sommerfeld);
            obj.scheme = scheme;
            obj.sommerfeld = sommerfeld;
        end        
        
    end
    
end

