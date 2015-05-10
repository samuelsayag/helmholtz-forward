classdef TestProblemSolver  < matlab.unittest.TestCase
    %TESTPROBLEMSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    methods  (Test)
        function test_build_matrix(testCase)
            param.m = 3;
            param.n = 3;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            solver = @(A, b) A\b;
            ps = ProblemSolver(param, scheme, solver);
            [ A, b, x ] = ps.solve();
            
%             full(A)
%             full(b)
%             full(x)            
            
            exp_A = [
                -3     1     0     1     0     0     0     0     0;
                 1    -3     1     0     1     0     0     0     0;
                 0     1    -3     0     0     1     0     0     0;
                 1     0     0    -3     1     0     1     0     0;
                 0     1     0     1    -3     1     0     1     0;
                 0     0     1     0     1    -3     0     0     1;
                 0     0     0     1     0     0    -3     1     0;
                 0     0     0     0     1     0     1    -3     1;
                 0     0     0     0     0     1     0     1    -3];            
            
            exp_b = [-10 -5 -10 -5 0 -5 -10 -5 -10]';
            
            exp_x = [
          26.666666666666647  34.999999999999972  26.666666666666647;
          34.999999999999964  46.666666666666615  34.999999999999964;
          26.666666666666647  34.999999999999972  26.666666666666643];
             
            t = {'vector b is not correct ', ...
                'Please recheck the class MatrixBuilder'};                        
            testCase.verifyEqual(full(b), exp_b, strjoin(t));            
            
            t = {'matrix A is not correct', ...
                'Please recheck the class MatrixBuilder'};            
            testCase.verifyEqual(full(A), exp_A, strjoin(t));            
                        
            t = {'vector x is not correct ', ...
                'Please recheck the class MatrixBuilder'};                        
            testCase.verifyEqual(full(x), exp_x, strjoin(t));                        
        end         
    end
    
end

