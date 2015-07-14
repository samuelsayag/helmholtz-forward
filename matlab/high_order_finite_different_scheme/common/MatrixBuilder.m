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
    end
    
    methods (Access = public)
        function obj = MatrixBuilder(basicScheme)
            obj = check_param(obj, basicScheme);
            obj.scheme = basicScheme;
            obj.mn = basicScheme.m * basicScheme.n;            
            % optimization zone
            p = gcp();
            obj.chunk_s = fix(obj.mn/p.NumWorkers);
            obj.stencil_size = 9;
        end
        
        function [A, b] = build(obj)
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
                last = [cartesian];
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
                [compl_idx, chunk] = fetchNext(f);
                for j = 1:size(chunk,1)
                    A( chunk{j,1} ) = chunk{j,2};
                    b( chunk{j,3} ) = chunk{j,4};                   
                end
            end                                    
            
% naive code without parallelisation for compare and debug
            
%             for i = 1:cartesian                
%                 [c_A, v_A, c_b, v_b] = obj.get_line(C1(i),C2(i));
%                 A(c_A) = v_A;
%                 b(c_b) = v_b;                                   
%             end                        
            
        end
    end
    
    methods (Access = private)
        
        function [chunk] = compute_chunk(obj, C1, C2, begi, endi)
            chunk = cell(endi - begi + 1, 4); cr = 1;
            for k = begi:endi                
                [c_A, v_A, c_b, v_b] = obj.get_line(C1(k),C2(k));           
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
        
        function [c_A, v_A, c_b, v_b] = get_line(obj, i, j)
            % function obj = get_line(obj, i, j)
            % i index of the lines ("vertically")
            % j index of the column ("horizontally")
            if i == 1
                if j == 1
                    [c_A, v_A, c_b, v_b] = obj.scheme.sw_pt( i, j );
                elseif j == obj.scheme.n
                    [c_A, v_A, c_b, v_b] = obj.scheme.se_pt( i, j );
                else
                    [c_A, v_A, c_b, v_b] = obj.scheme.s_pt( i, j );
                end
            elseif i == obj.scheme.m
               if j==1
                    [c_A, v_A, c_b, v_b] = obj.scheme.nw_pt( i, j );
               elseif j== obj.scheme.n
                    [c_A, v_A, c_b, v_b] = obj.scheme.ne_pt( i, j );
               else
                    [c_A, v_A, c_b, v_b] = obj.scheme.n_pt( i, j );
               end
            elseif j == 1
                [c_A, v_A, c_b, v_b] = obj.scheme.w_pt( i, j );
            elseif j == obj.scheme.n
                [c_A, v_A, c_b, v_b] = obj.scheme.e_pt( i, j );
            else
                [c_A, v_A, c_b, v_b] = obj.scheme.c_pt( i, j );
            end
            
        end
        
        function obj = check_param(obj, scheme)
            p = inputParser;

            schemes = {'BasicScheme'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));
            
            parse(p, scheme);            
        end 
        
    end   
end

