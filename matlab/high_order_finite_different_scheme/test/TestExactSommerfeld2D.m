classdef TestExactSommerfeld2D < matlab.unittest.TestCase                                
    %TESTORD6THHEMHOLTZ2D test the class ORD2NDHEMHOLTZ2D    
     
    methods (Test)
        
        function test_sx(testCase)
            h = 0.1;
            k = 10;
            theta = 0;
            somm = Ord2ndHelmholtz2D(k, h);
            scheme = ExactSommerfeld2D(h, k, theta, somm);
            t = {'not expected result for coefficient sx. '...
                'Please recheck the class ExactSommerfeld2D'};         
            
            coeff = 2 * 1i *(sin(h * k * cos(theta)));
            testCase.verifyEqual(scheme.sx, coeff, strjoin(t));
        end
        
        function test_sy(testCase)
            h = 0.1;
            k = 10;
            theta = pi/2;

            somm = Ord2ndHelmholtz2D(k, h);
            scheme = ExactSommerfeld2D(h, k, theta, somm);
            t = {'not expected result for coefficient sy. '...
                'Please recheck the class ExactSommerfeld2D'};         
            
            coeff = 2 * 1i *(sin(h * k * sin(theta)));
            testCase.verifyEqual(scheme.sy, coeff, strjoin(t));
        end
        
    end   
end

