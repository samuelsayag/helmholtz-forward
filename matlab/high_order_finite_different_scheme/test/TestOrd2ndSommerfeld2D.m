classdef TestOrd2ndSommerfeld2D < matlab.unittest.TestCase                                
    %TESTORD2NDSOMMERFELD2D test the class ORD6THSOMMERFELD2D    
     
    methods (Test)
        
        function test_sx(testCase)
            beta.x = 10;
            scheme = Ord2ndSommerfeld2D(0.1, beta);
            t = {'not expected result for coefficient sx. '...
                'Please recheck the class Ord2ndSommerfeld2D'};
            testCase.verifyEqual(scheme.sx, 2i, strjoin(t));
        end
        
        function test_sy(testCase)
            beta.y = 10;
            scheme = Ord2ndSommerfeld2D(0.1, beta);
            t = {'not expected result for coefficient sy. '...
                'Please recheck the class Ord2ndSommerfeld2D'};
            testCase.verifyEqual(scheme.sy, 2i, strjoin(t));
        end
        
    end   
end

