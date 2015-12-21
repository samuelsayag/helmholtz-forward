classdef ExactScheme2D  < CentralScheme
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
        % Optional parameters:
        % bessel_p1 = the first parameter it may have two meanings. If the
        % function choosen is bessel at exact theta this parameter is the
        % angle theta needed for computation (see: bessel_exact_theta 
        % function. If two parameter are given as it is the first bound for
        % the integral (see bessel_integral)
        % bessel_p2 = the second bound of the bessel integral.
        
            narginchk(2, 4)            
            obj = obj.check_param( k, h );
            
            if nargin == 2
                obj.bessel = obj.bessel_std;
            elseif nargin == 3
                obj = obj.check_theta( bessel_p1 );
                obj.bessel = @(x) obj.bessel_exact(x, bessel_p1);
            elseif nargin == 4
                obj = obj.check_integral( bessel_p1, bessel_p2 );
                obj.bessel = @(x) obj.bessel_integral(x, bessel_p1, bessel_p2);
            end            
        end        
        
        function a0 = a0(obj)
        % return A0 coefficient
            a0 = 4 * obj.bessel(obj.k * obj.h); 
        end
        
        function as = as(obj)
        % return the As coefficient
            obj.h; % dummy instruction    
            as = -1;
        end
        
        function ac = ac(obj)
        % return the Ac coefficient
            obj.h; % dummy instruction<
            ac = 0; 
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
        end
        
        function obj = check_integral(obj, bessel_p1, bessel_p2)
            p = inputParser;
            addRequired(p, 'bessel_p1', @isnumeric);
            addRequired(p, 'bessel_p2', @isnumeric);
            parse(p, bessel_p1, bessel_p2);                        
        end
        
    end 
    
end

