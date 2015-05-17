classdef TestDirichlet  < matlab.unittest.TestCase
    %TESTDIRICHLET set of test for the Dirichlet boundary problem
    
    methods (Test)
        
        function test_boundary(testCase)            
            % basic parameter of the simulation
            param.k = 5;
            param.h = 0.02;
            % definition of the area simulated
            param.a = 0; 
            param.b = 1;
            param.c = -1/2; 
            param.d = 1/2;
            param.m = (param.d - param.c)/param.h + 1;
            param.n = (param.b - param.a)/param.h + 1;

            % dirichlet function
            param.dirichlet = @(x,y) helm_sol1( x, y, param.k );            
            param.dirichlet = DirichletBuilder( param );
            
            param
            
            x = linspace(param.a-param.h, param.b+param.h, param.m + 2);
            y = linspace(param.d+param.h, param.c-param.h, param.n + 2);            
            [X,Y] = meshgrid(x,y);
            
            % build a matrix with just the helm func value on border and
            % zeros elswhere
            An = helm_sol1( X, Y, param.k );
            An(2:param.m+1, 2:param.n+1) = 0;
            
            % we build a similar matrix but withe the dirichlet builder
            % function the use the i,j index in the matrix to see if border
            % are calculated correctly while building the matrix.
            A = ones(size(An));
            A(2:param.m+1, 2:param.n+1) = 0;
            A = sparse(A);
            [row,col] = find(A);
            for i = 1:size(row,1)
                A(row(i),col(i)) = param.dirichlet( col(i)-1, (row(i)-1) );
            end
            
            nr = norm(An - full(A), Inf);
            testCase.assertTrue(nr < 1e-15, 'not equal');
        end
        
    end
    
end

