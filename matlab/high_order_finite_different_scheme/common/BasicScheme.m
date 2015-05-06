classdef BasicScheme
    %BASICSCHEME This class support the general functionality of a basic
    %9-point stencil scheme (which may be specialized in a 5-point
    %scheme).
    
    properties (SetAccess = public)
        % An instance of a scheme object that expose parameter of the
        % scheme in the form : A0 u_ij + As sig_s + Ac sig_c = 0
        % 2012 Erlangga, Turkel - ITERATIVE SCHEMES FOR HIGH ORDER COMPACT 
        % DISCRETIZATIONS
        % schemes = {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D',...
        % 'Ord6thHelmholtz2D'};
        scheme;
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
        param;
        % An alias short-cut that give more direct access to the dirichlet
        % function.
        dirichlet;
    end
    
    methods (Access = public)
        function obj = BasicScheme(param, scheme)
            narginchk(2, 2)
            obj = obj.check_param(param, scheme);            
            obj.param = param;
            obj.scheme = scheme;
            obj.dirichlet = param.dirichlet;
        end                
        
        function m = m(obj)
            m = obj.param.m;
        end
        
        function n = n(obj)
            n = obj.param.n;
        end
        
        function [c_A, v_A, c_b, v_b] = c_pt( obj, i, j )
            [c_A, c_b] = obj.c_pt_coordinate( i, j );
            [v_A, v_b] = obj.c_pt_value();
        end
        
        function [c_A, v_A, c_b, v_b] = n_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).
            [c_A, v_A, c_b, v_b] = obj.n_pt_dir( i, j );
        end
        
        function [c_A, v_A, c_b, v_b] = e_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).
            [c_A, v_A, c_b, v_b] = obj.e_pt_dir( i, j );
        end   
        
        function [c_A, v_A, c_b, v_b] = s_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).
            [c_A, v_A, c_b, v_b] = obj.s_pt_dir( i, j );
        end   
                
        function [c_A, v_A, c_b, v_b] = w_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer (Sommerfeld, Dirichlet,
            % Neumann, Robin...).
            [c_A, v_A, c_b, v_b] = obj.w_pt_dir( i, j );
        end
        
        function [c_A, v_A, c_b, v_b] = ne_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            [c_A, v_A, c_b, v_b] = obj.ne_pt_dir( i, j );        
        end
                
        function [c_A, v_A, c_b, v_b] = se_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            [c_A, v_A, c_b, v_b] = obj.se_pt_dir( i, j );                    
        end        
        
        function [c_A, v_A, c_b, v_b] = sw_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            [c_A, v_A, c_b, v_b] = obj.sw_pt_dir( i, j );                    
        end
        
        function [c_A, v_A, c_b, v_b] = nw_pt( obj, i, j )
            % TODO introduce a strategy pattern. The result return will be
            % that of a function pointer that may handle pure Dirichlet,
            % and variable combination with Sommerfeld + Dirichlet
            [c_A, v_A, c_b, v_b] = obj.nw_pt_dir( i, j );                    
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       PRIVATE PART OF THE CLASS                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Access = private)        
    
        function l = lin_lab( obj, i, j )
            % linear labelling matlab compliant:
            % http://www.mathworks.com/help/releases/R2014a/matlab/math/matrix-indexing.html
            d1 = obj.param.m * obj.param.n;
            l = (j-1) * d1 + i;
        end
        
        function l = label(obj, i,j)
            % considering that the area we solve the problem on is a square
            % of m * n points. This formula give the unique label of a
            % point so that the resulting matrix is k-diagonal.
            l = j + (obj.param.m - i) * obj.param.n;
        end
        
        function l = lab_n( obj, i, j )
            % give the label of a north point
            l = obj.label(i+1,j);
        end
                
        function l = lab_ne( obj, i, j )
            % give the label of a north east point
            l = obj.label(i+1,j+1);
        end        
                
        function l = lab_e( obj, i, j )
            % give the label of a north east point
            l = obj.label(i,j+1);
        end        
                
        function l = lab_se( obj, i, j )
            % give the label of a south east point
            l = obj.label(i-1,j+1);
        end
                
        function l = lab_s( obj, i, j )
            % give the label of a south point
            l = obj.label(i-1,j);
        end
                
        function l = lab_sw( obj, i, j )
            % give the label of a south west point
            l = obj.label(i-1,j-1);
        end
                
        function l = lab_w( obj, i, j )
            % give the label of a west point
            l = obj.label(i,j-1);
        end
                        
        function l = lab_nw( obj, i, j )
            % give the label of a north west point
            l = obj.label(i+1,j-1);
        end
        
        function [c_A, c_b] = c_pt_coordinate( obj, i, j )
            % provide the central stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(9,1); l = obj.label(i,j);                                    
            c_A(1) = obj.lin_lab(l, l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_n(i,j)); % north point
            c_A(3) = obj.lin_lab(l, obj.lab_ne(i,j)); % north east point
            c_A(4) = obj.lin_lab(l, obj.lab_e(i,j)); % east point
            c_A(5) = obj.lin_lab(l, obj.lab_se(i,j)); % south east point
            c_A(6) = obj.lin_lab(l, obj.lab_s(i,j)); % south point
            c_A(7) = obj.lin_lab(l, obj.lab_sw(i,j)); % south west point
            c_A(8) = obj.lin_lab(l, obj.lab_w(i,j)); % west point
            c_A(9) = obj.lin_lab(l, obj.lab_nw(i,j)); % north west point
            
            c_b = l; % vector b coordinate
        end 
        
        function [v_A, v_b] = c_pt_value( obj )
            % provide the central stencil points value
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b 
            v_A = zeros(9,1);            
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.ac; % north east point
            v_A(4) = obj.scheme.as; % east point
            v_A(5) = obj.scheme.ac; % south east point
            v_A(6) = obj.scheme.as; % south point
            v_A(7) = obj.scheme.ac; % south west point
            v_A(8) = obj.scheme.as; % west point
            v_A(9) = obj.scheme.ac; % north west point
            
            v_b = 0; % value is null in vector b for a central point
        end

        
        function l = dir_n( obj, i, j )
            % give the dirichlet value of a north point
            l = obj.dirichlet(i,j+1);
        end
                
        function l = dir_ne( obj, i, j )
            % give the dirichlet value of a north east point
            l = obj.dirichlet(i+1,j+1);
        end        
                
        function l = dir_e( obj, i, j )
            % give the dirichlet value of a north east point
            l = obj.dirichlet(i+1,j);
        end        
                
        function l = dir_se( obj, i, j )
            % give the dirichlet value of a south east point
            l = obj.dirichlet(i+1,j-1);
        end
                
        function l = dir_s( obj, i, j )
            % give the dirichlet value of a south point
            l = obj.dirichlet(i,j-1);
        end
                
        function l = dir_sw( obj, i, j )
            % give the dirichlet value of a south west point
            l = obj.dirichlet(i-1,j-1);
        end
                
        function l = dir_w( obj, i, j )
            % give the dirichlet value of a west point
            l = obj.dirichlet(i-1,j);
        end
                        
        function l = dir_nw( obj, i, j )
            % give the dirichlet value of a north west point
            l = obj.dirichlet(i-1,j+1);
        end        
        
        function [c_A, v_A, c_b, v_b] = n_pt_dir( obj, i, j )
            [c_A, c_b] = obj.n_pt_coordinate( i, j );
            [v_A, v_b] = obj.n_pt_value_dirichlet( i, j );
        end        
        
        function [c_A, c_b] = n_pt_coordinate( obj, i, j )
            % provide the north side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1); l = obj.label(i,j);
            c_A(1) = obj.lin_lab(l,l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_e(i,j)); % east point
            c_A(3) = obj.lin_lab(l, obj.lab_se(i,j)); % south east point
            c_A(4) = obj.lin_lab(l, obj.lab_s(i,j)); % south point
            c_A(5) = obj.lin_lab(l, obj.lab_sw(i,j)); % south west point
            c_A(6) = obj.lin_lab(l, obj.lab_w(i,j)); % west point
            
            c_b = l; % vector b coordinate            
        end
        
        function [v_A, v_b] = n_pt_value_dirichlet( obj, i, j )
            % provide the north side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % east point
            v_A(3) = obj.scheme.ac; % south east point
            v_A(4) = obj.scheme.as; % south point
            v_A(5) = obj.scheme.ac; % south west point
            v_A(6) = obj.scheme.as; % west point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.bc * obj.dir_nw(i,j) ...
                + obj.scheme.bs * obj.dir_n(i,j) ...
                + obj.scheme.bc * obj.dir_ne(i,j) );
        end        
        
        function [c_A, v_A, c_b, v_b] = e_pt_dir( obj, i, j )
            [c_A, c_b] = obj.e_pt_coordinate( i, j );
            [v_A, v_b] = obj.e_pt_value_dirichlet( i, j );
        end        
        
        function [c_A, c_b] = e_pt_coordinate( obj, i, j )
            % provide the east side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1); l = obj.label(i,j);
            c_A(1) = obj.lin_lab(l,l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_n(i,j)); % north point
            c_A(3) = obj.lin_lab(l, obj.lab_s(i,j)); % south point
            c_A(4) = obj.lin_lab(l, obj.lab_sw(i,j)); % south west point
            c_A(5) = obj.lin_lab(l, obj.lab_w(i,j)); % west point
            c_A(6) = obj.lin_lab(l, obj.lab_nw(i,j)); % north west point
            
            c_b = l; % vector b coordinate            
        end
        
        function [v_A, v_b] = e_pt_value_dirichlet( obj, i, j )
            % provide the east side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.as; % south point
            v_A(4) = obj.scheme.ac; % south west point
            v_A(5) = obj.scheme.as; % west point
            v_A(6) = obj.scheme.ac; % north west point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.bc * obj.dir_ne(i,j)...
                + obj.scheme.bs * obj.dir_e(i,j)...
                + obj.scheme.bc * obj.dir_se(i,j) );
        end
    
        function [c_A, v_A, c_b, v_b] = s_pt_dir( obj, i, j )
            [c_A, c_b] = obj.s_pt_coordinate( i, j );
            [v_A, v_b] = obj.s_pt_value_dirichlet( i, j );
        end        
        
        function [c_A, c_b] = s_pt_coordinate( obj, i, j )
            % provide the south side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1); l = obj.label(i,j);
            c_A(1) = obj.lin_lab(l,l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_n(i,j)); % north point
            c_A(3) = obj.lin_lab(l, obj.lab_ne(i,j)); % north east point
            c_A(4) = obj.lin_lab(l, obj.lab_e(i,j)); % east point
            c_A(5) = obj.lin_lab(l, obj.lab_w(i,j)); % west point
            c_A(6) = obj.lin_lab(l, obj.lab_nw(i,j)); % north west point
            
            c_b = l; % vector b coordinate            
        end
        
        function [v_A, v_b] = s_pt_value_dirichlet( obj, i, j )
            % provide the south side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.ac; % north east point
            v_A(4) = obj.scheme.as; % east point
            v_A(5) = obj.scheme.as; % west point
            v_A(6) = obj.scheme.ac; % north west point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.bc * obj.dir_se(i,j)...
                + obj.scheme.bs * obj.dir_s(i,j)...
                + obj.scheme.bc * obj.dir_sw(i,j) );
        end
            
        function [c_A, v_A, c_b, v_b] = w_pt_dir( obj, i, j )
            [c_A, c_b] = obj.w_pt_coordinate( i, j );
            [v_A, v_b] = obj.w_pt_value_dirichlet( i, j );
        end        
        
        function [c_A, c_b] = w_pt_coordinate( obj, i, j )
            % provide the south side stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b            
            c_A = zeros(6,1); l = obj.label(i,j);
            c_A(1) = obj.lin_lab(l,l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_n(i,j)); % north point
            c_A(3) = obj.lin_lab(l, obj.lab_ne(i,j)); % north east point
            c_A(4) = obj.lin_lab(l, obj.lab_e(i,j)); % east point
            c_A(5) = obj.lin_lab(l, obj.lab_se(i,j)); % south east point
            c_A(6) = obj.lin_lab(l, obj.lab_s(i,j)); % south point            
            c_b = l; % vector b coordinate            
        end
        
        function [v_A, v_b] = w_pt_value_dirichlet( obj, i, j )
            % provide the west side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(6,1);
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.ac; % north east point
            v_A(4) = obj.scheme.as; % east point
            v_A(5) = obj.scheme.ac; % south east point
            v_A(6) = obj.scheme.as; % south point
            
            % !!! the sign here is strongly linked to the instance of the
            % scheme and the way it is written !!!
            v_b = - ( obj.scheme.bc * obj.dir_sw(i,j)...
                + obj.scheme.bs * obj.dir_w(i,j)...
                + obj.scheme.bc * obj.dir_nw(i,j) );
        end
                    
        function [c_A, v_A, c_b, v_b] = ne_pt_dir( obj, i, j )
            [c_A, c_b] = obj.ne_pt_coordinate( i, j );
            [v_A, v_b] = obj.ne_pt_value_dirichlet( i, j );
        end        
            
        function [c_A, v_A, c_b, v_b] = se_pt_dir( obj, i, j )
            [c_A, c_b] = obj.se_pt_coordinate( i, j );
            [v_A, v_b] = obj.se_pt_value_dirichlet( i, j );
        end        
                    
        function [c_A, v_A, c_b, v_b] = sw_pt_dir( obj, i, j )
            [c_A, c_b] = obj.sw_pt_coordinate( i, j );
            [v_A, v_b] = obj.sw_pt_value_dirichlet( i, j );
        end        
            
        function [c_A, v_A, c_b, v_b] = nw_pt_dir( obj, i, j )
            [c_A, c_b] = obj.nw_pt_coordinate( i, j );
            [v_A, v_b] = obj.nw_pt_value_dirichlet( i, j );
        end        
        
        function [c_A, c_b] = ne_pt_coordinate( obj, i, j )
            % provide the NE corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1); l = obj.label(i,j);                                    
            c_A(1) = obj.lin_lab(l, l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_s(i,j)); % south point
            c_A(3) = obj.lin_lab(l, obj.lab_sw(i,j)); % south west point
            c_A(4) = obj.lin_lab(l, obj.lab_w(i,j)); % west point
            
            c_b = l; % vector b coordinate
            
        end
        
        function [c_A, c_b] = se_pt_coordinate( obj, i, j )
            % provide the SE corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1); l = obj.label(i,j);                                    
            c_A(1) = obj.lin_lab(l, l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_n(i,j)); % north point
            c_A(3) = obj.lin_lab(l, obj.lab_w(i,j)); % west point
            c_A(4) = obj.lin_lab(l, obj.lab_nw(i,j)); % north west point
            
            c_b = l; % vector b coordinate            
        end
        
        function [c_A, c_b] = sw_pt_coordinate( obj, i, j )
            % provide the SW corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1); l = obj.label(i,j);                                    
            c_A(1) = obj.lin_lab(l, l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_n(i,j)); % north point
            c_A(3) = obj.lin_lab(l, obj.lab_ne(i,j)); % north east point
            c_A(4) = obj.lin_lab(l, obj.lab_e(i,j)); % east point
            
            c_b = l; % vector b coordinate            
        end
        
        function [c_A, c_b] = nw_pt_coordinate( obj, i, j )
            % provide the NW corner stencil points coordinate as matlab 
            % compliant linear labelling in the matrix A and vector b for
            % the equation Ax=b
            % return:
            %   c_A: the linear coordinate of the different stencil points
            %   in the matrix A
            %   c_b: the linear coordinate in the vector b
            c_A = zeros(4,1); l = obj.label(i,j);                                    
            c_A(1) = obj.lin_lab(l, l); % central point            
            c_A(2) = obj.lin_lab(l, obj.lab_e(i,j)); % east point
            c_A(3) = obj.lin_lab(l, obj.lab_se(i,j)); % south east point
            c_A(4) = obj.lin_lab(l, obj.lab_s(i,j)); % south point
            
            c_b = l; % vector b coordinate            
        end
        
        function [v_A, v_b] = ne_pt_value_dirichlet( obj, i, j )
            % provide the NE side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % south point
            v_A(3) = obj.scheme.ac; % south west point
            v_A(4) = obj.scheme.as; % west point
            
            % value is dirichlet for other point 
            v_b = -(obj.scheme.bc * obj.dir_nw(i,j)...
                + obj.scheme.bs * obj.dir_n(i,j) ...
                + obj.scheme.bc * obj.dir_ne(i,j)....
                + obj.scheme.bs * obj.dir_e(i,j)....
                + obj.scheme.bc * obj.dir_se(i,j) ); 
        end
        
        function [v_A, v_b] = se_pt_value_dirichlet( obj, i, j )
            % provide the NE side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b        
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.as; % west point
            v_A(4) = obj.scheme.ac; % north west point

            % value is dirichlet for other point 
            v_b = -(obj.scheme.bc * obj.dir_ne(i,j)...
                + obj.scheme.bs * obj.dir_e(i,j)...
                + obj.scheme.bc * obj.dir_se(i,j) ...
                + obj.scheme.bs * obj.dir_s(i,j)...
                + obj.scheme.bc * obj.dir_sw(i,j) );             
        end
        
        function [v_A, v_b] = sw_pt_value_dirichlet( obj, i, j )
            % provide the SW side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b        
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % north point
            v_A(3) = obj.scheme.ac; % north east point
            v_A(4) = obj.scheme.as; % east point

            % value is dirichlet for other point 
            v_b = -( obj.scheme.bc * obj.dir_nw(i,j)...
                + obj.scheme.bs * obj.dir_n(i,j)...
                + obj.scheme.bc * obj.dir_ne(i,j)...
                + obj.scheme.bs * obj.dir_e(i,j)...
                + obj.scheme.bc * obj.dir_se(i,j) );            
        end
        
        function [v_A, v_b] = nw_pt_value_dirichlet( obj, i, j )
            % provide the NW side stencil points value + dirichlet
            % return:
            %   v_A: the coeff value of the different stencil points
            %   in the matrix A
            %   v_b: the coeff value in the vector b        
            v_A = zeros(4,1);            
            v_A(1) = obj.scheme.a0; % central point            
            v_A(2) = obj.scheme.as; % east point
            v_A(3) = obj.scheme.ac; % south east point
            v_A(4) = obj.scheme.as; % south point

            % value is dirichlet for other point 
            v_b = -( obj.scheme.bc * obj.dir_sw(i,j)...
                + obj.scheme.bs * obj.dir_w(i,j)...
                + obj.scheme.bc * obj.dir_nw(i,j)...
                + obj.scheme.bs * obj.dir_n(i,j)...
                + obj.scheme.bc * obj.dir_ne(i,j) );             
        end
        
        function obj = check_param(obj, param, scheme)
            p = inputParser;

            schemes = {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D',...
                'Ord6thHelmholtz2D'};
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
      
    end    
end