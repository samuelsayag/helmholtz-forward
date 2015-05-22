classdef TestExactHemholtz2D < matlab.unittest.TestCase                                
    %TESTORD2NDHEMHOLTZ2D test the class ORD2NDHEMHOLTZ2D    
    
    methods (Test)
        
        function test_a0_coeff(testCase)
            scheme = ExactScheme2D(10, 0.1);
            t = {'not expected result for coefficient a0. '...
                'Please recheck the class ExactScheme2D'};
            e = scheme.a0 - (-3.060790746231866);
            testCase.verifyLessThan(e, 1e-15, strjoin(t));            
        end
                
        function test_as_coeff(testCase)
            scheme = ExactScheme2D(10, 0.1);
            t = {'not expected result for coefficient as. '...
                'Please recheck the class ExactScheme2D'};
            testCase.verifyEqual(scheme.as, 1, strjoin(t));
        end
                
        function test_ac_coeff(testCase)
            scheme = ExactScheme2D(10, 0.1);
            t = {'not expected result for coefficient ac. '...
                'Please recheck the class ExactScheme2D'};
            testCase.verifyEqual(scheme.ac, 0, strjoin(t));
        end
        
    end
    
end

