classdef MatrixBuilder
    %MATRIXBUILDER this class take some parameter of a simulation to
    %instanciate and is responsible for building the matrix and vector b
    %such that : Ax=b with x the solution of the problem.
    
    properties (SetAccess = private)
        scheme; % the basic scheme object
        
        % space and time optimization zone
        % for sparse matrix preallocation of non zeros elements
        stencil_size; 
        % size of data packet processed by each threads
        chunk_s;
    end
    
    properties (SetAccess = private)
        A; % the matrix that is built
        b; % the vector b that is built
        mn;
        
        % the building type of the matrix ('sequential', 'parall_1' 
        % 'parall_2')
        parall = [];        
    end
    
    methods (Access = public)
        function obj = MatrixBuilder(basicScheme, parall)
            obj = check_param(obj, basicScheme);
            obj.scheme = basicScheme;
            obj.mn = basicScheme.m * basicScheme.n; 

            % define the parallelisation mode for building the matrix
            obj.parall = parall;
            % optimization zone
            obj.stencil_size = 9;            
        end
        
        function [A, b] = build(obj)
%==========================================================================
%    NO PARALELISATION METHOD (maintain for debug and performance study)
%==========================================================================                                  
          % without parallelisation
          if strcmp(obj.parall, 'sequential')             
            [A, b] = naive_building(obj);
          elseif any(strcmp(obj.parall, {'parall_1', 'parall_2'}))
            p = gcp();
            obj.chunk_s = fix(obj.mn/p.NumWorkers);
%             obj.stencil_size = 9;               
            if strcmp(obj.parall,'parall_1')
                % type 1 parallisation worker/chunk a
                [A, b] = obj.parallel_1();               
            end
            if strcmp(obj.parall,'parall_2')
%                 type two parallisation parfor/variable slicing
                [A, b] = obj.parallel_2();
            end            
          else
            error('Parallelisation mode: sequential, parall_1, parall_2');
          end
%==========================================================================
%    DECOMMENT THIS CODE WATHEVER PARALLELISATION METHOD CHOOSEN
%==========================================================================            
%             p = gcp();
%             obj.chunk_s = fix(obj.mn/p.NumWorkers);
%             obj.stencil_size = 9;
%==========================================================================
%             FIRST PARALELISATION METHOD
%==========================================================================            
%             % type 1 parallisation worker/chunk
%             [A, b] = obj.parallel_1();
%==========================================================================
%             SECOND PARALELISATION METHOD
%==========================================================================            
            % type two parallisation parfor/variable slicing
%             [A, b] = obj.parallel_2();
        end
    end
    
    methods (Access = private)
        
        function [A, b] = parallel_2(obj)
            % Try to optimize time of construction of the matrix A by
            % pre-allocating space for the non zeros numbers
            nzmax = obj.scheme.m * obj.stencil_size;
            A = spalloc(obj.mn, obj.mn, nzmax);
            b = spalloc(obj.mn,1, 2 * obj.mn);

            % Create two vectors of i and j indexes over the modelled grid  
            [C1, C2] = meshgrid(...
                linspace(1,obj.scheme.m, obj.scheme.m),...
                linspace(obj.scheme.n, 1, obj.scheme.n));            
            C1 = C1(:); C2 = C2(:);

            % create a matrix (3 * mn) matrix. Each line is a 3 value 
            % vector (final line label in A and B, i, j)
            % The sort is required by the parfor loop !!
            cartesian = size(C1,1); 
            line_labels = LabelHandler( obj.scheme.m, obj.scheme.n, C1, C2 );
            tab_label = sortrows([line_labels, C1, C2], 1);

            % Avoid unnecessary copy of tab_label in the parfor by providing 
            % just a function pointer. A broadcast bariable 
            look_up = @(l,c) tab_label(l,c); 

            parfor l = 1:cartesian       
                % get the coordinate of the value in the line (c_A) and the
                % values (v_A) the value in the vector b (v_b)                
                [c_A, v_A, ~, v_b] = obj.get_line('ll', look_up(l,2), look_up(l,3));
                % fill the matrix and the vector at line corresponding to 
                % the computed label                
                A(l,:) = sparse(1, c_A, v_A, 1, obj.mn);
                b(l,:) = v_b;
            end
            
        end
        
        function [A, b] = parallel_1(obj)
            % try to optimize time of construction of the matrix A by
            % preallocating space for the non zeros numbers
            nzmax = obj.scheme.m * obj.stencil_size;
            A = spalloc(obj.mn, obj.mn, nzmax);
            b = spalloc(obj.mn,1, 2 * obj.mn);
            
            % create two vectors of i and j indexes over the modeled grid
            [C1, C2] = obj.get_index_vectors();                                    
            cartesian = size(C1,1); 
            last = [];
            if (mod(cartesian, obj.chunk_s))
                last = cartesian;
            end
            idx_chunk = [obj.chunk_s : obj.chunk_s : cartesian last];
            
            % compute the lines for each chunk
            func = @(w,x,y,z) obj.compute_chunk(w,x,y,z);
            begi = 1; N = size(idx_chunk,2); p = gcp();            
            for idx = 1:N 
                endi = idx_chunk(idx);
                f(idx) = parfeval(p, func, 1, C1, C2, begi, endi);
                begi = begi + obj.chunk_s;
            end            
            % bring out the results and fill the matrix
            for i = 1:N                
                [~, chunk] = fetchNext(f);
                for j = 1:size(chunk,1)
                    A( chunk{j,1} ) = chunk{j,2};
                    b( chunk{j,3} ) = chunk{j,4};                   
                end
            end                                    
            
        end
        
        function [A, b] = naive_building(obj)
            nzmax = obj.scheme.m * obj.stencil_size;
            A = spalloc(obj.mn, obj.mn, nzmax);
            b = spalloc(obj.mn,1, 2 * obj.mn);
            
            % create two vectors of i and j indexes over the modeled grid
            [C1, C2] = obj.get_index_vectors();                                    
            cartesian = size(C1,1); 
            
            % naive code without parallelisation for compare and debug            
            for i = 1:cartesian                
                [c_A, v_A, c_b, v_b] = obj.get_line('ml', C1(i),C2(i));
                A(c_A) = v_A;
                b(c_b) = v_b;                                   
            end                        
        
        end
        
        function [chunk] = compute_chunk(obj, C1, C2, begi, endi)
            chunk = cell(endi - begi + 1, 4); cr = 1;
            for k = begi:endi                
                [c_A, v_A, c_b, v_b] = obj.get_line('ml', C1(k),C2(k));           
                chunk(cr,1:4) = {c_A, v_A, c_b, v_b};
                cr = cr + 1;
            end                        
        end
        
        function [C1, C2] = get_index_vectors(obj)
            % build indexing for a one loop for (instead of a nested loop)
            [C1, C2] = meshgrid(...
                linspace(1,obj.scheme.m, obj.scheme.m),...
                linspace(obj.scheme.n, 1, obj.scheme.n));
            C1 = C1(:); C2 = C2(:);        
        end        
        
        function [c_A, v_A, c_b, v_b] = get_line(obj, type_coord, i, j)
            if i == 1
                if j == 1
                    [c_A, v_A, c_b, v_b] = obj.scheme.sw_pt( type_coord, i, j );
                elseif j == obj.scheme.n
                    [c_A, v_A, c_b, v_b] = obj.scheme.se_pt( type_coord, i, j );
                else
                    [c_A, v_A, c_b, v_b] = obj.scheme.s_pt( type_coord, i, j );
                end
            elseif i == obj.scheme.m
               if j==1
                    [c_A, v_A, c_b, v_b] = obj.scheme.nw_pt( type_coord, i, j );
               elseif j== obj.scheme.n
                    [c_A, v_A, c_b, v_b] = obj.scheme.ne_pt( type_coord, i, j );
               else
                    [c_A, v_A, c_b, v_b] = obj.scheme.n_pt( type_coord, i, j );
               end
            elseif j == 1
                [c_A, v_A, c_b, v_b] = obj.scheme.w_pt( type_coord, i, j );
            elseif j == obj.scheme.n
                [c_A, v_A, c_b, v_b] = obj.scheme.e_pt( type_coord, i, j );
            else
                [c_A, v_A, c_b, v_b] = obj.scheme.c_pt( type_coord, i, j );
            end
            
        end
        
        function obj = check_param(obj, scheme)
            p = inputParser;

            schemes = {'BasicScheme', 'BasicScheme2'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));
            
            parse(p, scheme);            
        end          
    end   
end

