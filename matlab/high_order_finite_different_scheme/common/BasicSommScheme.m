classdef BasicSommScheme
    %BASICSOMMSCHEME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
                obj.scheme.as * obj.sommerfeld.sy; % central point            
            v_A(2) = obj.scheme.as - obj.scheme.ac * obj.sommerfeld.sy; % east point
            v_A(3) = 2 * obj.scheme.ac; % south east point
            v_A(4) = 2 * obj.scheme.as; % south point
            v_A(5) = 2 * obj.scheme.ac; % south west point
            v_A(6) = obj.scheme.as - obj.scheme.ac * obj.sommerfeld.sy; % west point            
        end
        
        function v_A = e_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 - ...
                obj.scheme.as * obj.sommerfeld.sx; % central point
            v_A(2) = obj.scheme.as - obj.scheme.ac * obj.sommerfeld.sx; % north point
            v_A(3) = obj.scheme.as - obj.scheme.ac * obj.sommerfeld.sx; % south point
            v_A(4) = 2 * obj.scheme.ac; % south west point
            v_A(5) = 2 * obj.scheme.as; % west point            
            v_A(6) = 2 * obj.scheme.ac; % north west point
        end
        
        function v_A = s_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 + ...
                obj.scheme.as * obj.sommerfeld.sy; % central point
            v_A(2) = 2 * obj.scheme.as; % north point
            v_A(3) = 2 * obj.scheme.ac; % north east point
            v_A(4) = obj.scheme.as + obj.scheme.ac * obj.sommerfeld.sy; % east point
            v_A(5) = obj.scheme.as + obj.scheme.ac * obj.sommerfeld.sy; % west point
            v_A(6) = 2 * obj.scheme.ac; % north west point
        end        
        
        function v_A = w_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 + ...
                obj.scheme.as * obj.sommerfeld.sx; % central point
            v_A(2) = obj.scheme.as + obj.scheme.ac * obj.sommerfeld.sx; % north point
            v_A(3) = 2 * obj.scheme.ac; % north east point
            v_A(4) = 2 * obj.scheme.as; % east point
            v_A(5) = 2 * obj.scheme.ac; % south east point
            v_A(6) = obj.scheme.as + obj.scheme.ac * obj.sommerfeld.sx; % south point
        end        
        
        function v_A = ne_pt_som_dir( obj )
            v_A = []; 
        end
        
        function v_A = se_pt_som_dir( obj )
            v_A = [];
        end
        
        function v_A = sw_pt_som_dir( obj )
            v_A = [];
        end        
        
        function v_A = nw_pt_som_dir( obj )
            v_A = [];
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

