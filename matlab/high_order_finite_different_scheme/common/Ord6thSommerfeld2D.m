classdef Ord6thSommerfeld2D
    %ORD6THHELMHOLTZSOMMERFELD2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        h;
        beta;
        scheme;
    end
    
    methods 
        function obj = Ord6thSommerfeld2D( h, beta, scheme )
        % Ord6thSommerfeld2D
        % h: the the step of the grid
        % beta: parameter of the formula (pi - k^2 for instance)
            narginchk(3, 3)
            obj = obj.check_param( h, beta, scheme );            
        end
        
        function v_A = n_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 - obj.scheme.as * obj.s0; % central point            
            v_A(2) = obj.scheme.as; % east point
            v_A(3) = obj.scheme.ac; % south east point
            v_A(4) = 2 * obj.scheme.as; % south point
            v_A(5) = obj.scheme.ac; % south west point
            v_A(6) = obj.scheme.as; % west point            
        end
        
        function v_A = e_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 - obj.scheme.as * obj.s0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.as; % south point
            v_A(4) = obj.scheme.ac; % south west point
            v_A(5) = 2 * obj.scheme.as; % west point
            v_A(6) = obj.scheme.ac; % north west point
        end
        
        function v_A = s_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 + obj.scheme.as * obj.s0; % central point            
            v_A(2) = 2 * obj.scheme.as; % north point
            v_A(3) = obj.scheme.ac; % north east point
            v_A(4) = obj.scheme.as; % east point
            v_A(5) = obj.scheme.as; % west point
            v_A(6) = obj.scheme.ac; % north west point
        end        
        
        function v_A = w_pt( obj )
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0 + obj.scheme.as * obj.s0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.ac; % north east point
            v_A(4) = 2 * obj.scheme.as; % east point
            v_A(5) = obj.scheme.ac; % south east point
            v_A(6) = obj.scheme.as; % south point
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
    
    methods (Access = private)
        
        function out = out(obj)
            obj.h; % dummy instruction
            out = 1;
        end
        
        function in = in(obj)
            obj.h; % dummy instruction
            in = 1;            
        end        
        
        function ctr = s0(obj)
            obj.h; % dummy instruction
            ctr = 1;            
        end        
        
        function obj = check_param(obj, h, beta, scheme)                                 
            p = inputParser;           
            
            addRequired(p, 'h', @isnumeric);   
            addRequired(p, 'beta', @isnumeric);     
            schemes = {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D',...
                'Ord6thHelmholtz2D'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));            
            
            parse( p, h, beta, scheme );            

            obj.k = k;
            obj.h = h;
            obj.scheme = scheme;
        end        
    end
    
end

