classdef ExactScheme2D
    %EXACTSCHEME2D This class contain the code for the so called exacte
    %scheme the appear in two of the following articles.
    %
    %Reference:
    %- A NEW FINITE DIFFERENCE METHOD FOR THE HELMHOLTZ EQUATION USING 
    % SYMBOLIC COMPUTATION, LARRY A. LAMBE, RICHARD LUCZAK, 
    % AND JOHN W. NEHRBASS, 2003.
    %- EXACT FINITE DIFFERENCE SCHEMES FOR SOLVING HELMHOLTZ EQUATION 
    % AT ANY WAVENUMBER, YAU SHU WONG AND GUANGRUI LI, 2011. 
   
    properties (SetAccess = private)
        bessel_std = @(x) besselj( 0, x );
        bessel_integral = @(x, a, b) bessel_integral( x, a, b );
        bessel_exact = @(x, theta) bessel_exact_theta( x, theta );
    end
    
    properties (SetAccess = public)
        k = [];
        h = [];
        bessel = [];
    end    
    

    methods (Access = public)    
        
        function obj = ExactScheme2D(k, h, bessel_p1, bessel_p2)
        % Ord2ndHelmholtz
        % k = the k of an Helmholtz equation 
        % h = the step length (basic division of the grid)
        % theta = (optional) necessary to compute the exact theta function
        % instead the the bessel integral.
            narginchk(2, 3)
            if nargin == 2
                obj.bessel = obj.bessel_std;
            elseif nargin == 3
                obj = obj.check_theta(obj, bessel_p1);
            elseif nargin == 4
                obj = obj.check_param(k, h);
            end
        end        
        
    end
    
    methods (Access = private)
        
        function obj = check_param(obj, k, h)
            p = inputParser;
            addRequired(p, 'k', @isnumeric);
            addRequired(p, 'h', @isnumeric);        
            parse(p, k, h);
            obj.k = k;
            obj.h = h;            
        end
        
        function obj = check_theta(obj, theta)
            p = inputParser;
            addRequired(p, 'theta', @isnumeric);
            parse(p, theta);
            obj.bessel = @(x) obj.bessel_exact(x, theta);
        end
        
    end 
    
end

