classdef BasicScheme
    %BASICSCHEME This class support the general functionality of a basic
    %9-point stencil scheme (which may be specialized in a 5-point
    %scheme).
    
    properties (SetAccess = private)
       dir_dir = {'dirichlet', 'dirichlet'}
    end
    
    properties (SetAccess = private)
        % An instance of a scheme object that expose parameter of the
        % scheme in the form : A0 u_ij + As sig_s + Ac sig_c = 0
        % 2012 Erlangga, Turkel - ITERATIVE SCHEMES FOR HIGH ORDER COMPACT 
        % DISCRETIZATIONS
        % schemes = {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D',...
        % 'Ord6thHelmholtz2D'};
        scheme = [];
        % Sommerfeld constraint under the form of a scheme with coefficient
        % just as for the central scheme.
        sommerfeld  = [];
        % A struct that fields are as follow:
        % param.h: length of the step in the unit chosen ex: 0.002 = 2mm if
        %   the meter is chosen as a unit.
        % param.k: the wave number k = 2 * pi * f / c (f = frequency,
        %   c = celerity of the wave)
        % param.m: the number of point along the y axis (number of line)
        % param.n: the number of point along the x axis (number of column)
        % param.dirichlet: the dirichlet function as pointer of function
        %   that give a value by passing it the line and column index. Ex:
        %   dirichlet = @(i,j) some_function....
        % param.north = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the north side of the area.
        % param.east = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the east side of the area.
        % param.south = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the south side of the area.
        % param.west = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the west side of the area.
        % 
        % Note: it will not be possible to choose Sommerfeld type
        % boundaries along all the side. At least one must be of Dirichlet
        % type.
        param  = [];
        % An alias short-cut that give more direct access to the dirichlet
        % function.
        dirichlet  = [];
    end
    
    methods (Access = public)
        function obj = BasicScheme(param, scheme, sommerfeld)
            narginchk(2, 3);
            obj = obj.check_param(param, scheme);            
            obj.param = param;
            obj.scheme = scheme;
            obj.dirichlet = param.dirichlet;
            if nargin == 3
                obj = obj.check_sommerfeld(sommerfeld);
            end
            obj = obj.check_boundary();
        end                
        
        function m = m(obj)
            m = obj.param.m;
        end
        
        function n = n(obj)
            n = obj.param.n;
        end
        
        function [c_A, v_A, c_b, v_b] = c_pt( obj, type_coord, i, j )
            l = obj.label(i,j);
            f_coord = obj.get_coordinate_func(type_coord, l);
            c_A = obj.c_pt_coordinate( f_coord, i, j );
            c_b = l;
            [v_A, v_b] = obj.c_pt_value();
        end
        
        function [c_A, v_A, c_b, v_b] = n_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).
            
            if strcmp(obj.param.north, 'sommerfeld')
                [c_A, v_A, c_b, v_b] = obj.n_pt_som( type_coord, i, j );                
            else
                [c_A, v_A, c_b, v_b] = obj.n_pt_dir( type_coord, i, j );
            end
        end
        
        function [c_A, v_A, c_b, v_b] = e_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).            
            if strcmp(obj.param.east, 'sommerfeld')
                [c_A, v_A, c_b, v_b] = obj.e_pt_som( type_coord, i, j );
            else
                [c_A, v_A, c_b, v_b] = obj.e_pt_dir( type_coord, i, j );
            end
        end   
        
        function [c_A, v_A, c_b, v_b] = s_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).
            if strcmp(obj.param.south, 'sommerfeld')
                [c_A, v_A, c_b, v_b] = obj.s_pt_som( type_coord, i, j );
            else
                [c_A, v_A, c_b, v_b] = obj.s_pt_dir( type_coord, i, j );
            end
        end   
                
        function [c_A, v_A, c_b, v_b] = w_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).
            if strcmp(obj.param.west, 'sommerfeld')
                [c_A, v_A, c_b, v_b] = obj.w_pt_som( type_coord, i, j );
            else
                [c_A, v_A, c_b, v_b] = obj.w_pt_dir( type_coord, i, j );
            end
        end
        
        function [c_A, v_A, c_b, v_b] = ne_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            t = strcmp({obj.param.north, obj.param.east}, ...
                obj.dir_dir);
            if all(t) % all dirichlet
                [c_A, v_A, c_b, v_b] = obj.ne_pt_dir( type_coord, i, j );            
            elseif ~any(t) % all sommerfeld
                [c_A, v_A, c_b, v_b] = obj.ne_pt_som_som( type_coord, i, j );
            elseif strcmp(obj.param.north, 'sommerfeld') % som_dir
                [c_A, v_A, c_b, v_b] = obj.ne_pt_som_dir( type_coord, i, j );
            else % dir_som
                [c_A, v_A, c_b, v_b] = obj.ne_pt_dir_som( type_coord, i, j );
            end        
        end
                
        function [c_A, v_A, c_b, v_b] = se_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            t = strcmp({obj.param.south, obj.param.east}, ...
                obj.dir_dir);
            if all(t) % all dirichlet
                [c_A, v_A, c_b, v_b] = obj.se_pt_dir( type_coord, i, j );                    
            elseif ~any(t) % all sommerfeld
                [c_A, v_A, c_b, v_b] = obj.se_pt_som_som( type_coord, i, j );                
            elseif strcmp(obj.param.south, 'sommerfeld') % som_dir
                [c_A, v_A, c_b, v_b] = obj.se_pt_som_dir( type_coord, i, j );
            else % dir_som
                [c_A, v_A, c_b, v_b] = obj.se_pt_dir_som( type_coord, i, j );
            end            
        end        
        
        function [c_A, v_A, c_b, v_b] = sw_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            t = strcmp({obj.param.south, obj.param.west}, ...
                obj.dir_dir);
            if all(t) % all dirichlet
                [c_A, v_A, c_b, v_b] = obj.sw_pt_dir( type_coord, i, j );                                    
            elseif ~any(t) % all sommerfeld
                [c_A, v_A, c_b, v_b] = obj.sw_pt_som_som( type_coord, i, j );                
            elseif strcmp(obj.param.south, 'sommerfeld') % som_dir
                [c_A, v_A, c_b, v_b] = obj.sw_pt_som_dir( type_coord, i, j );
            else % dir_som
                [c_A, v_A, c_b, v_b] = obj.sw_pt_dir_som( type_coord, i, j );
            end                        
        end
        
        function [c_A, v_A, c_b, v_b] = nw_pt( obj, type_coord, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            t = strcmp({obj.param.north, obj.param.west}, ...
                obj.dir_dir);
            if all(t) % all dirichlet
                [c_A, v_A, c_b, v_b] = obj.nw_pt_dir( type_coord, i, j );                                    
            elseif ~any(t) % all sommerfeld
                [c_A, v_A, c_b, v_b] = obj.nw_pt_som_som( type_coord, i, j );                
            elseif strcmp(obj.param.north, 'sommerfeld') % som_dir
                [c_A, v_A, c_b, v_b] = obj.nw_pt_som_dir( type_coord, i, j );
            else % dir_som
                [c_A, v_A, c_b, v_b] = obj.nw_pt_dir_som( type_coord, i, j );
            end                                    
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       PRIVATE PART OF THE CLASS                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    methods(Static, Access = private)
        function c_A = c_pt_coordinate( f_coord, i, j )
            % provide the central stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            
            c_A = zeros(9,1);                                              
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i+1,j); % north point
            c_A(3) = f_coord(i+1,j+1); % north east point
            c_A(4) = f_coord(i,j+1); % east point
            c_A(5) = f_coord(i-1,j+1); % south east point
            c_A(6) = f_coord(i-1,j); % south point
            c_A(7) = f_coord(i-1,j-1); % south west point
            c_A(8) = f_coord(i,j-1); % west point
            c_A(9) = f_coord(i+1,j-1); % north west point
        end
        
        function c_A = n_pt_coordinate( f_coord, i, j )
            % provide the north side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1); 
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i,j+1); % east point
            c_A(3) = f_coord(i-1,j+1); % south east point
            c_A(4) = f_coord(i-1,j); % south point
            c_A(5) = f_coord(i-1,j-1); % south west point
            c_A(6) = f_coord(i,j-1); % west point                                  
        end
                       
        function c_A = e_pt_coordinate( f_coord, i, j )
            % provide the east side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1);             
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i+1,j); % north point
            c_A(3) = f_coord(i-1,j); % south point
            c_A(4) = f_coord(i-1,j-1); % south west point
            c_A(5) = f_coord(i,j-1); % west point
            c_A(6) = f_coord(i+1,j-1); % north west point            
        end                
        
        function c_A = s_pt_coordinate( f_coord, i, j )
            % provide the south side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1);
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i+1,j); % north point
            c_A(3) = f_coord(i+1,j+1); % north east point
            c_A(4) = f_coord(i,j+1); % east point
            c_A(5) = f_coord(i,j-1); % west point
            c_A(6) = f_coord(i+1,j-1); % north west point            
        end                
        
        function c_A = w_pt_coordinate( f_coord, i, j )
            % provide the south side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1);            
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i+1,j); % north point
            c_A(3) = f_coord(i+1,j+1); % north east point
            c_A(4) = f_coord(i,j+1); % east point
            c_A(5) = f_coord(i-1,j+1); % south east point
            c_A(6) = f_coord(i-1,j); % south point           
        end 
        function c_A = ne_pt_coordinate( f_coord, i, j )
            % provide the NE corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1);             
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i-1,j); % south point
            c_A(3) = f_coord(i-1,j-1); % south west point
            c_A(4) = f_coord(i,j-1); % west point
        end
        
        function c_A = se_pt_coordinate( f_coord, i, j )
            % provide the SE corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1);            
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i+1,j); % north point
            c_A(3) = f_coord(i,j-1); % west point
            c_A(4) = f_coord(i+1,j-1); % north west point            
        end
        
        function c_A = sw_pt_coordinate( f_coord, i, j )
            % provide the SW corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1);             
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i+1,j); % north point
            c_A(3) = f_coord(i+1,j+1); % north east point
            c_A(4) = f_coord(i,j+1); % east point           
        end
        
        function c_A = nw_pt_coordinate( f_coord, i, j )
            % provide the NW corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1);
            c_A(1) = f_coord(i,j); % central point    
            c_A(2) = f_coord(i,j+1); % east point
            c_A(3) = f_coord(i-1,j+1); % south east point
            c_A(4) = f_coord(i-1,j); % south point
        end        
    end
    
    methods (Access = private)        
    
        function l = lin_lab( obj, i, j )
            % linear labelling matlab compliant:
            % http://www.mathworks.com/help/releases/R2014a/matlab/math/matrix-indexing.html
            d1 = obj.param.m * obj.param.n;
            l = (j-1) * d1 + i;
        end
        
        function l = label(obj, i, j)
            % considering that the area we solve the problem on is a square
            % of m * n points. This formula give the unique label of a
            % point so that the resulting matrix is k-diagonal.
            l = j + (obj.param.m - i) * obj.param.n;
        end
        
        function coord_func = get_coordinate_func(obj, type_coord, line_label)
            % obtain two type of function depending of the type of
            % coordinate we want our results in:
            % matrix_linear is one coordinate to describe a cell in a
            % matrix
            % line_label is ones coordinate to descibe a cell in a line of
            % a matrix
            if strcmp(type_coord, 'ml')
                coord_func = @(i, j) ...
                    obj.lin_lab(line_label, obj.label(i,j));
                return 
            elseif strcmp(type_coord, 'll')
                coord_func = @(i, j, line_label) obj.label(i,j);
                return
            else
                error('Unknown function type of coordinate');
            end
        end

        function [v_A, v_b] = c_pt_value( obj )
            % provide the central stencil points value
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b

            v_A = zeros(9,1);            
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.n; % north point
            v_A(3) = obj.scheme.ne; % north east point
            v_A(4) = obj.scheme.e; % east point
            v_A(5) = obj.scheme.se; % south east point
            v_A(6) = obj.scheme.s; % south point
            v_A(7) = obj.scheme.sw; % south west point
            v_A(8) = obj.scheme.w; % west point
            v_A(9) = obj.scheme.nw; % north west point

            v_b = 0; % value is null in vector b for a central point
        end

        function d = dir_n( obj, i, j )
            % give the dirichlet value of a north point
            d = obj.dirichlet(j,i+1);
        end
                
        function d = dir_ne( obj, i, j )
            % give the dirichlet value of a north east point
            d = obj.dirichlet(j+1,i+1);
        end        
                
        function d = dir_e( obj, i, j )
            % give the dirichlet value of a north east point
            d = obj.dirichlet(j+1,i);
        end        
                
        function d = dir_se( obj, i, j )
            % give the dirichlet value of a south east point
            d = obj.dirichlet(j+1,i-1);
        end
                
        function d = dir_s( obj, i, j )
            % give the dirichlet value of a south point
            d = obj.dirichlet(j,i-1);
        end
                
        function d = dir_sw( obj, i, j )
            % give the dirichlet value of a south west point
            d = obj.dirichlet(j-1,i-1);
        end
                
        function d = dir_w( obj, i, j )
            % give the dirichlet value of a west point
            d = obj.dirichlet(j-1,i);
        end
                        
        function d = dir_nw( obj, i, j )
            % give the dirichlet value of a north west point
            d = obj.dirichlet(j-1,i+1);
        end        
        
        
        function l = dirichlet_coordinate(obj, i, j )
            l = obj.label(i,j);
        end       
       
        
        function [c_A, v_A, c_b, v_b] = n_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.n_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.n_pt_value_dirichlet( i, j );
        end         
                
        function [c_A, v_A, c_b, v_b] = e_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.e_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.e_pt_value_dirichlet( i, j );
        end           
    
        function [c_A, v_A, c_b, v_b] = s_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.s_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.s_pt_value_dirichlet( i, j );
        end        
            
        function [c_A, v_A, c_b, v_b] = w_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.w_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.w_pt_value_dirichlet( i, j );
        end                
        
        function [v_A, v_b] = n_pt_value_dirichlet( obj, i, j )
            % provide the north side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.e; % east point
            v_A(3) = obj.scheme.se; % south east point
            v_A(4) = obj.scheme.s; % south point
            v_A(5) = obj.scheme.sw; % south west point
            v_A(6) = obj.scheme.w; % west point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.nw * obj.dir_nw(i,j) ...
                + obj.scheme.n * obj.dir_n(i,j) ...
                + obj.scheme.ne * obj.dir_ne(i,j) );
        end             
        
        function [v_A, v_b] = e_pt_value_dirichlet( obj, i, j )
            % provide the east side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.n; % north point
            v_A(3) = obj.scheme.s; % south point
            v_A(4) = obj.scheme.sw; % south west point
            v_A(5) = obj.scheme.w; % west point
            v_A(6) = obj.scheme.nw; % north west point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.ne * obj.dir_ne(i,j)...
                + obj.scheme.e * obj.dir_e(i,j)...
                + obj.scheme.se * obj.dir_se(i,j) );
        end
        
        function [v_A, v_b] = s_pt_value_dirichlet( obj, i, j )
            % provide the south side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.n; % north point
            v_A(3) = obj.scheme.ne; % north east point
            v_A(4) = obj.scheme.e; % east point
            v_A(5) = obj.scheme.w; % west point
            v_A(6) = obj.scheme.nw; % north west point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.se * obj.dir_se(i,j)...
                + obj.scheme.s * obj.dir_s(i,j)...
                + obj.scheme.sw * obj.dir_sw(i,j) );
        end
        
        function [v_A, v_b] = w_pt_value_dirichlet( obj, i, j )
            % provide the west side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.n; % north point
            v_A(3) = obj.scheme.ne; % north east point
            v_A(4) = obj.scheme.e; % east point
            v_A(5) = obj.scheme.se; % south east point
            v_A(6) = obj.scheme.s; % south point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.sw * obj.dir_sw(i,j)...
                + obj.scheme.w * obj.dir_w(i,j)...
                + obj.scheme.nw * obj.dir_nw(i,j) );
        end
                    
        function [c_A, v_A, c_b, v_b] = ne_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.ne_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.ne_pt_value_dirichlet( i, j );
        end        
            
        function [c_A, v_A, c_b, v_b] = se_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.se_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.se_pt_value_dirichlet( i, j );
        end        
                    
        function [c_A, v_A, c_b, v_b] = sw_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.sw_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.sw_pt_value_dirichlet( i, j );
        end        
            
        function [c_A, v_A, c_b, v_b] = nw_pt_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.nw_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.nw_pt_value_dirichlet( i, j );
        end        
        

        
        function [v_A, v_b] = ne_pt_value_dirichlet( obj, i, j )
            % provide the NE side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.s; % south point
            v_A(3) = obj.scheme.sw; % south west point
            v_A(4) = obj.scheme.w; % west point
            
            % value is dirichlet for other point 
            v_b = -(obj.scheme.nw * obj.dir_nw(i,j)...
                + obj.scheme.n * obj.dir_n(i,j) ...
                + obj.scheme.ne * obj.dir_ne(i,j)....
                + obj.scheme.e * obj.dir_e(i,j)....
                + obj.scheme.se * obj.dir_se(i,j) ); 
        end        
        
        function [v_A, v_b] = se_pt_value_dirichlet( obj, i, j )
            % provide the NE side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b        
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.n; % north point
            v_A(3) = obj.scheme.w; % west point
            v_A(4) = obj.scheme.nw; % north west point

            % value is dirichlet for other point 
            v_b = -(obj.scheme.ne * obj.dir_ne(i,j)...
                + obj.scheme.e * obj.dir_e(i,j)...
                + obj.scheme.se * obj.dir_se(i,j) ...
                + obj.scheme.s * obj.dir_s(i,j)...
                + obj.scheme.sw * obj.dir_sw(i,j) );             
        end
        
        function [v_A, v_b] = sw_pt_value_dirichlet( obj, i, j )
            % provide the SW side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b        
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.n; % north point
            v_A(3) = obj.scheme.ne; % north east point
            v_A(4) = obj.scheme.e; % east point

            % value is dirichlet for other point 
            v_b = -( obj.scheme.nw * obj.dir_nw(i,j)...
                + obj.scheme.w * obj.dir_w(i,j)...
                + obj.scheme.sw * obj.dir_sw(i,j)...
                + obj.scheme.s * obj.dir_s(i,j)...
                + obj.scheme.se * obj.dir_se(i,j) );            
        end
        
        function [v_A, v_b] = nw_pt_value_dirichlet( obj, i, j )
            % provide the NW side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b        
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.c; % central point            
            v_A(2) = obj.scheme.e; % east point
            v_A(3) = obj.scheme.se; % south east point
            v_A(4) = obj.scheme.s; % south point

            % value is dirichlet for other point 
            v_b = -( obj.scheme.sw * obj.dir_sw(i,j)...
                + obj.scheme.w * obj.dir_w(i,j)...
                + obj.scheme.nw * obj.dir_nw(i,j)...
                + obj.scheme.n * obj.dir_n(i,j)...
                + obj.scheme.ne * obj.dir_ne(i,j) );             
        end                

        function [v_b] = n_half_ne_corner_dirichlet( obj, i, j )
            % The north half of the north east corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.nw * obj.dir_nw(i,j)...
                + obj.scheme.n * obj.dir_n(i,j) ...
                + obj.scheme.ne * obj.dir_ne(i,j) ); 
        end        

        function [v_b] = e_half_ne_corner_dirichlet( obj, i, j )
            % The east half of the north east corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.ne() * obj.dir_ne(i,j)....
                + obj.scheme.e() * obj.dir_e(i,j)....
                + obj.scheme.se() * obj.dir_se(i,j) ); 
        end        
        
        function [v_b] = e_half_se_corner_dirichlet( obj, i, j )
            % The east half of the south east corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b        
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.ne() * obj.dir_ne(i,j)...
                + obj.scheme.e() * obj.dir_e(i,j)...
                + obj.scheme.se() * obj.dir_se(i,j) );             
        end
        
        function [v_b] = s_half_se_corner_dirichlet( obj, i, j )
            % The south half of the south east corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b        
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.se() * obj.dir_se(i,j) ...
                + obj.scheme.s() * obj.dir_s(i,j)...
                + obj.scheme.sw() * obj.dir_sw(i,j) );             
        end
        
        function [v_b] = s_half_sw_corner_dirichlet( obj, i, j )
            % The south half of the south west corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b        
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.sw() * obj.dir_sw(i,j)...
                + obj.scheme.s() * obj.dir_s(i,j)...
                + obj.scheme.se() * obj.dir_se(i,j) );            
        end
        
        function [v_b] = w_half_sw_corner_dirichlet( obj, i, j )
            % The west half of the south west corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b        
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.nw() * obj.dir_nw(i,j)...
                + obj.scheme.w() * obj.dir_w(i,j)...
                + obj.scheme.sw() * obj.dir_sw(i,j) );            
        end
        
        function [v_b] = w_half_nw_corner_dirichlet( obj, i, j )
            % The west half of the north west corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b        
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.sw() * obj.dir_sw(i,j)...
                + obj.scheme.w() * obj.dir_w(i,j)...
                + obj.scheme.nw() * obj.dir_nw(i,j) );             
        end                
        
        function [v_b] = n_half_nw_corner_dirichlet( obj, i, j )
            % The west half of the north west corner is provided here as a
            % dirichlet boundary.
            % return:
            %   v_b: the coeff value in the vector b        
            
            % value is dirichlet for other point 
            v_b = -( obj.scheme.nw() * obj.dir_nw(i,j)...
                + obj.scheme.n() * obj.dir_n(i,j)...
                + obj.scheme.ne() * obj.dir_ne(i,j) );             
        end                        
        
        function [c_A, v_A, c_b, v_b] = n_pt_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.n_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.n_pt_value_som();
        end         
                
        function [c_A, v_A, c_b, v_b] = e_pt_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.e_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.e_pt_value_som();
        end           
    
        function [c_A, v_A, c_b, v_b] = s_pt_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.s_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.s_pt_value_som();
        end        
            
        function [c_A, v_A, c_b, v_b] = w_pt_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.w_pt_coordinate( f_coord, i, j );
            c_b = obj.dirichlet_coordinate( i, j );
            [v_A, v_b] = obj.w_pt_value_som();
        end         
        
        function [v_A, v_b] = n_pt_value_som( obj )
            % provide the north side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = obj.sommerfeld.n_pt();            
            % No dirichlet the b vector receive 0
            v_b = 0;
        end        
        
        function [v_A, v_b] = e_pt_value_som( obj )
            % provide the east side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = obj.sommerfeld.e_pt();                        
            % No dirichlet the b vector receive 0
            v_b = 0;
        end
        
        function [v_A, v_b] = s_pt_value_som( obj )
            % provide the south side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = obj.sommerfeld.s_pt();                                    
            % No dirichlet the b vector receive 0
            v_b = 0;
        end
        
        function [v_A, v_b] = w_pt_value_som( obj )
            % provide the west side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = obj.sommerfeld.w_pt();                                    
            % No dirichlet the b vector receive 0
            v_b = 0;
        end        
                        
        function [c_A, v_A, c_b, v_b] = ne_pt_som_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.ne_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.n_half_ne_pt();
            v_b = obj.e_half_ne_corner_dirichlet( i, j );            
        end
        
        function [c_A, v_A, c_b, v_b] = ne_pt_dir_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.ne_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.e_half_ne_pt();
            v_b = obj.n_half_ne_corner_dirichlet( i, j );
        end        
                
        function [c_A, v_A, c_b, v_b] = se_pt_som_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.se_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.s_half_se_pt();
            v_b = obj.e_half_se_corner_dirichlet( i, j );
        end           
    
        function [c_A, v_A, c_b, v_b] = se_pt_dir_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.se_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.e_half_se_pt();
            v_b = obj.s_half_se_corner_dirichlet( i, j );
        end        
            
        function [c_A, v_A, c_b, v_b] = sw_pt_som_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.sw_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.s_half_sw_pt();
            v_b = obj.w_half_sw_corner_dirichlet( i, j );
        end         
            
        function [c_A, v_A, c_b, v_b] = sw_pt_dir_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.sw_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.w_half_sw_pt();
            v_b = obj.s_half_sw_corner_dirichlet( i, j );
        end         
            
        function [c_A, v_A, c_b, v_b] = nw_pt_som_dir( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.nw_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.n_half_nw_pt();
            v_b = obj.w_half_nw_corner_dirichlet( i, j );
        end         
            
        function [c_A, v_A, c_b, v_b] = nw_pt_dir_som( obj, type_coord, i, j )
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.nw_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.w_half_nw_pt();
            v_b = obj.n_half_nw_corner_dirichlet( i, j );
        end                                 
                
        function [c_A, v_A, c_b, v_b] = ne_pt_som_som( obj, type_coord, i, j )
            % provide the north east corner sommerfeld coefficient
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.ne_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.ne_pt();            
            % No dirichlet the b vector receive 0
            v_b = 0;              
        end        
                        
        function [c_A, v_A, c_b, v_b] = se_pt_som_som( obj, type_coord, i, j )
            % provide the south east corner sommerfeld coefficient
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.se_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.se_pt();            
            % No dirichlet the b vector receive 0
            v_b = 0;              
        end
                
        function [c_A, v_A, c_b, v_b] = sw_pt_som_som( obj, type_coord, i, j )
            % provide the south west corner sommerfeld coefficient
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.sw_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.sw_pt();            
            % No dirichlet the b vector receive 0
            v_b = 0;              
        end
                
        function [c_A, v_A, c_b, v_b] = nw_pt_som_som( obj, type_coord, i, j )
            % provide the north west corner sommerfeld coefficient
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            f_coord = obj.get_coordinate_func(type_coord, obj.label(i,j));
            c_A = obj.nw_pt_coordinate( f_coord, i, j );            
            c_b = obj.dirichlet_coordinate( i, j );            
            v_A = obj.sommerfeld.nw_pt();            
            % No dirichlet the b vector receive 0
            v_b = 0;              
        end
        
        function obj = check_param(obj, param, scheme)
            p = inputParser;

%             schemes = {'Ord2ndHelmholtz2D',... 
%                 'Ord4thHelmholtz2D', 'Ord4thHelmholtz2D_2', ...
%                 'Ord6thHelmholtz2D', 'Ord6thHelmholtz2D_2', ...
%                 'ExactScheme2D', 'Poisson2D'};
            schemes = {'NinePtStencil'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));            
            
            function res = valide_param(x)
                validateattributes( x, {'struct'}, {'nonempty'});
                res = isfield(x, 'dirichlet');
                validateattributes( x.dirichlet, {'function_handle'}, {'nonempty'});
                res = isfield(x, 'm') && res;
                validateattributes( x.m, {'double'}, {'nonempty', 'integer'});                
                res = isfield(x, 'n') && res;
                validateattributes( x.n, {'double'}, {'nonempty','integer'});
            end
            addRequired(p, 'param', @(x)valide_param(x));                 
            
            parse(p, scheme, param);                        
        end
      
        function obj = check_sommerfeld(obj, sommerfeld)
            p = inputParser;
            schemes = {'BasicSommScheme', 'BasicSommScheme2'};
            addRequired(p, 'sommerfeld', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));            
            parse(p, sommerfeld);
            obj.sommerfeld = sommerfeld;
        end
        
        function obj = check_boundary( obj )
            if ~ isfield(obj.param, 'north')
                obj.param.north = 'dirichlet';
            end
            if ~ isfield(obj.param, 'east')
                obj.param.east = 'dirichlet';
            end
            if ~ isfield(obj.param, 'south')
                obj.param.south = 'dirichlet';
            end
            if ~ isfield(obj.param, 'west')
                obj.param.west = 'dirichlet';
            end
            
            valid_bounds = {'dirichlet','sommerfeld'};
            validatestring( obj.param.north, valid_bounds );
            validatestring( obj.param.east, valid_bounds );
            validatestring( obj.param.south, valid_bounds );
            validatestring( obj.param.west, valid_bounds );
            
            cond = strcmp(obj.param.north, 'sommerfeld') ...
                && strcmp(obj.param.south, 'sommerfeld') ...
                && strcmp(obj.param.east, 'sommerfeld') ...
                && strcmp(obj.param.west, 'sommerfeld');
            if cond
                error('Not all side may be sommerfeld.');
            end

            cond = (strcmp(obj.param.north, 'sommerfeld') ...
                || strcmp(obj.param.south, 'sommerfeld') ...
                || strcmp(obj.param.east, 'sommerfeld') ...
                || strcmp(obj.param.west, 'sommerfeld')) ...
            && isempty(obj.sommerfeld);            
            if cond 
                error('A sommerfeld scheme must be set.');
            end            
        end
        
    end    
end