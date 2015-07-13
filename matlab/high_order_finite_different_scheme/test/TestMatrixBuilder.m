classdef TestMatrixBuilder < matlab.unittest.TestCase                                
    %TESTMATRIXBUILDER test the class ORD2NDHEMHOLTZ2D    
    
    methods (Test)
        
        function test_build_matrix(testCase)
            param.m = 3;
            param.n = 3;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            mb = MatrixBuilder(bs);
            [A, b] = mb.build();
            
%             full(A)
%             full(b)'
            
            exp_A = -[...
            -3     1     0     1     0     0     0     0     0;
             1    -3     1     0     1     0     0     0     0; 
             0     1    -3     0     0     1     0     0     0;
             1     0     0    -3     1     0     1     0     0;
             0     1     0     1    -3     1     0     1     0;
             0     0     1     0     1    -3     0     0     1;
             0     0     0     1     0     0    -3     1     0;
             0     0     0     0     1     0     1    -3     1;
             0     0     0     0     0     1     0     1    -3];
         
            exp_b = -[-10 -5 -10 -5  0 -5 -10 -5 -10]';
            
            t = {'vector b is not correct ', ...
                'Please recheck the class MatrixBuilder'};                        
            testCase.verifyEqual(full(b), exp_b, strjoin(t));            
            
            t = {'matrix A is not correct', ...
                'Please recheck the class MatrixBuilder'};            
            testCase.verifyEqual(full(A), exp_A, strjoin(t));            
        end                      
    end    
end

