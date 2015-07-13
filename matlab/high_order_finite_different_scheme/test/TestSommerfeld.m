classdef TestSommerfeld  < matlab.unittest.TestCase
    %TESTSOMMERFELD tests of the different Sommerfeld scheme
    
    methods (Test)
        function test_sommerfeld_exact_scheme_33(TestCase)
            clear variables; close all; clc;

            % modeled solution
            theor = @(x, y, k, theta) helm_sol2_2D( k, theta, x, y);
            theta = pi/4;

            % basic parameter of the simulation
            param.k = 5;
            param.h = 0.5;
            % definition of the area we simulate in it
            param.a = 0; 
            param.b = 1;
            param.c = 0; 
            param.d = 1;
            param.m = (param.d - param.c)/param.h + 1;
            param.n = (param.b - param.a)/param.h + 1;

            % dirichlet function
            param.dirichlet = @(x,y) theor( x, y, param.k , theta);
            scheme = ExactScheme2D(param.k, param.h);
            param.east = 'sommerfeld';
            beta = param.k;
            sommerfeld = ExactSommerfeld2D( param.h, beta, theta, scheme);

            % define the solver
            solver = @(A, b) A\b;

            param
            ps = ProblemSolver(param, scheme, solver, sommerfeld);
            [ A, b, sol ] = ps.solve();            
        end
    end
    
end

