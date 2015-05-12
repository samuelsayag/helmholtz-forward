classdef TestOrd6thHemholtz2D < matlab.unittest.TestCase                                
    %TESTORD6THHEMHOLTZ2D test the class ORD2NDHEMHOLTZ2D    
     
    methods (Test)
        
        function test_n_pt(testCase)
            scheme = Ord6thHelmholtz2D(sqrt(30), 1, 3);
            t = {'not expected result for coefficient a0. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.a0, 57/3, strjoin(t));
        end
        
    end   
end

