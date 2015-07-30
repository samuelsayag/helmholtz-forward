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
            v_A(1) = obj.scheme.c - obj.scheme.n .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.e - obj.scheme.ne .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % east point
            v_A(3) = obj.scheme.se - obj.scheme.ne .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south east point
            v_A(4) = obj.scheme.s -  obj.scheme.n .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south point
            v_A(5) = obj.scheme.sw - obj.scheme.nw .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south west point
            v_A(6) = obj.scheme.w - obj.scheme.nw .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % west point            
        end
        
        function v_A = e_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.c - obj.scheme.e .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % central point
            v_A(2) = obj.scheme.n - obj.scheme.ne .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % north point
            v_A(3) = obj.scheme.s - obj.scheme.se .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % south point
            v_A(4) = obj.scheme.sw - obj.scheme.se .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south west point
            v_A(5) = obj.scheme.w - obj.scheme.e .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % west point            
            v_A(6) = obj.scheme.nw - obj.scheme.ne .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % north west point
        end
        
        function v_A = s_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.c - obj.scheme.s .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.n - obj.scheme.s .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north point
            v_A(3) = obj.scheme.ne - obj.scheme.se .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north east point
            v_A(4) = obj.scheme.e - obj.scheme.se .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % east point
            v_A(5) = obj.scheme.w - obj.scheme.sw .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % west point
            v_A(6) = obj.scheme.nw - obj.scheme.sw .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north west point
        end        
        
        function v_A = w_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.c - obj.scheme.w .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.n - obj.scheme.nw .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % north point
            v_A(3) = obj.scheme.ne - obj.scheme.nw .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north east point
            v_A(4) = obj.scheme.e - obj.scheme.w .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % east point
            v_A(5) = obj.scheme.se - obj.scheme.sw .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % south east point
            v_A(6) = obj.scheme.s - obj.scheme.sw .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % south point
        end        
        
        function v_A = n_half_ne_pt( obj )
            % Sommerfeld is assumed on the North side but not on the East
            % side (east side is Dirichlet).
            v_A = zeros(4,1);                                    
            v_A(1) = obj.scheme.c - obj.scheme.n .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.s -  obj.scheme.n .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south point
            v_A(3) = obj.scheme.sw - obj.scheme.nw .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south west point
            v_A(4) = obj.scheme.w - obj.scheme.nw .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % west point            
        end
        
        function v_A = e_half_ne_pt( obj )  
            % Sommerfeld is assumed on the East side but not on the North
            % side (Dirichlet).            
            v_A = zeros(4,1);              
            v_A(1) = obj.scheme.c - obj.scheme.e .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % central point
            v_A(2) = obj.scheme.s - obj.scheme.se .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % south point
            v_A(3) = obj.scheme.sw - obj.scheme.se .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south west point
            v_A(4) = obj.scheme.w - obj.scheme.e .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % west point                   
        end
        
        function v_A = e_half_se_pt( obj )
            % Sommerfeld is assumed on the East side but not on the South
            % side (Dirichlet).              
            v_A = zeros(4,1);
            v_A(1) = obj.scheme.c - obj.scheme.e .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % central point
            v_A(2) = obj.scheme.n - obj.scheme.ne .* obj.sommerfeld.sox ./ obj.sommerfeld.sf; % north point
            v_A(3) = obj.scheme.w - obj.scheme.e .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % west point            
            v_A(4) = obj.scheme.nw - obj.scheme.ne .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % north west point                   
        end
        
        function v_A = s_half_se_pt( obj )   
            % Sommerfeld is assumed on the South side but not on the East
            % side (Dirichlet).            
            v_A = zeros(4,1);          
            v_A(1) = obj.scheme.c - obj.scheme.s .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.n - obj.scheme.s .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north point
            v_A(3) = obj.scheme.w - obj.scheme.sw .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % west point
            v_A(4) = obj.scheme.nw - obj.scheme.sw .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north west point            
        end
        
        function v_A = s_half_sw_pt( obj )
            % Sommerfeld is assumed on the South side but not on the West
            % side (Dirichlet).                        
            v_A = zeros(4,1);
            v_A(1) = obj.scheme.c - obj.scheme.s .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.n - obj.scheme.s .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north point
            v_A(3) = obj.scheme.ne - obj.scheme.se .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north east point
            v_A(4) = obj.scheme.e - obj.scheme.se .* obj.sommerfeld.soy ./ obj.sommerfeld.sb; % east point
        end        
        
        function v_A = w_half_sw_pt( obj )            
            v_A = zeros(4,1);  
            v_A(1) = obj.scheme.c - obj.scheme.w .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.n - obj.scheme.nw .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % north point
            v_A(3) = obj.scheme.ne - obj.scheme.nw .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % north east point
            v_A(4) = obj.scheme.e - obj.scheme.w .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % east point
        end        
        
        function v_A = w_half_nw_pt( obj )            
            v_A = zeros(4,1);  
            v_A(1) = obj.scheme.c - obj.scheme.w .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % central point
            v_A(2) = obj.scheme.e - obj.scheme.w .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % east point
            v_A(3) = obj.scheme.se - obj.scheme.sw .* obj.sommerfeld.sf ./ obj.sommerfeld.sb; % south east point
            v_A(4) = obj.scheme.s - obj.scheme.sw .* obj.sommerfeld.sox ./ obj.sommerfeld.sb; % south point        
        end         
        
        function v_A = n_half_nw_pt( obj ) 
            % Sommerfeld is assumed on the North side but not on the West
            % side (east side is Dirichlet).            
            v_A = zeros(4,1);
            v_A(1) = obj.scheme.c - obj.scheme.n .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % central point            
            v_A(2) = obj.scheme.e - obj.scheme.ne .* obj.sommerfeld.soy ./ obj.sommerfeld.sf; % east point
            v_A(3) = obj.scheme.se - obj.scheme.ne .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south east point
            v_A(4) = obj.scheme.s -  obj.scheme.n .* obj.sommerfeld.sb ./ obj.sommerfeld.sf; % south point
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

%             schemes = {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D',...
%                 'Ord6thHelmholtz2D', 'ExactScheme2D'};
%             addRequired(p, 'scheme', ...
%                 @(x)validateattributes( x, schemes, {'nonempty'}));            
            schemes = {'NinePtStencil'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));
            
%             sommerfelds = {'Ord6thSommerfeld2D', 'Ord2ndSommerfeld2D', ...
%                 'ExactSommerfeld2D'};
%             addRequired(p, 'sommerfeld', ...
%                 @(x)validateattributes( x, sommerfelds, {'nonempty'}));
            schemes = {'SommerfeldScheme'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));                        
            
            parse(p, scheme, sommerfeld);
            obj.scheme = scheme;
            obj.sommerfeld = sommerfeld;
        end        
        
    end
    
end

