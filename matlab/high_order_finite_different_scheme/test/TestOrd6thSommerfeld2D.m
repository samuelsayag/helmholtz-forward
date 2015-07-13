classdef TestOrd6thSommerfeld2D < matlab.unittest.TestCase                                
    %TESTORD6THSOMMERFELD2D test the class ORD6THSOMMERFELD2D    
     
    methods (Test)
        
        
        function test_sx(testCase)
            beta.x = 10;
            scheme = Ord6thSommerfeld2D(0.1, beta);
            t = {'not expected result for coefficient sx. '...
                'Please recheck the class Ord6thSommerfeld2D'};
            coeff = 2 * (1-1/6+1/120) * 1i;
            testCase.verifyEqual(scheme.sox, coeff, strjoin(t));
        end
        
        function test_sy(testCase)
            beta.y = 10;
            scheme = Ord6thSommerfeld2D(0.1, beta);
            t = {'not expected result for coefficient sy. '...
                'Please recheck the class Ord6thSommerfeld2D'};
            coeff = 2 * (1-1/6+1/120) * 1i;
            testCase.verifyEqual(scheme.soy, coeff, strjoin(t));
        end
        
    end   
end

