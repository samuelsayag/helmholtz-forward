classdef MatrixBuilder
    %MATRIXBUILDER this class take some parameter of a simulation to
    %instanciate and is responsible for building the matrix and vector b
    %such that : Ax=b with x the solution of the problem.
    
    properties (SetAccess = public)
        scheme; % the basic scheme object
    end
    
    properties (SetAccess = private)
        A; % the matrix that is built
        b; % the vector b that is built
    end
    
    methods (Access = public)
        function obj = MatrixBuilder(basicScheme)
            obj.scheme = basicScheme;
            mn = basicScheme.m * basicScheme.n;
            obj.A = sparse(mn,mn);
            obj.b = sparse(mn,1);
        end
        
        function [A, b] = build(obj)
            for i = 1:obj.scheme.m
                for j = 1:obj.scheme.n
                    obj = obj.get_line(i,j);
                end
            end
            A = obj.A;
            b = obj.b;
        end
    end
    
    methods (Access = private)
        
        function obj = get_line(obj, i, j)
            if i == 1
                if j==1
                    [c_A, v_A, c_b, v_b] = obj.scheme.sw_pt( i, j );
                elseif j==n
                    [c_A, v_A, c_b, v_b] = obj.schemes.nw_pt( i, j );
                else
                    [c_A, v_A, c_b, v_b] = obj.schemes.w_pt( i, j );
                end
            elseif i == m
               if j==1
                    [c_A, v_A, c_b, v_b] = obj.scheme.se_pt( i, j );
               elseif j==n
                    [c_A, v_A, c_b, v_b] = obj.scheme.ne_pt( i, j );
               else
                    [c_A, v_A, c_b, v_b] = obj.scheme.e_pt( i, j );
               end
            elseif j == 1
                [c_A, v_A, c_b, v_b] = obj.scheme.s_pt( i, j );
            elseif j == m
                [c_A, v_A, c_b, v_b] = obj.scheme.n_pt( i, j );
            else
                [c_A, v_A, c_b, v_b] = obj.scheme.c_pt( i, j );
            end
            
            obj.A(c_A) = v_A;
            obj.A(c_b) = v_b;
        end
    end
    
end

