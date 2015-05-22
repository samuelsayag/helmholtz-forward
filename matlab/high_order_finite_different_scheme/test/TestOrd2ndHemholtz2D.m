classdef TestOrd2ndHemholtz2D < matlab.unittest.TestCase                                
    %TESTORD2NDHEMHOLTZ2D test the class ORD2NDHEMHOLTZ2D    
    properties
        h;
        k;
    end
 
    methods(TestMethodSetup)
        function create_h_k(testCase)
            testCase.h = 0.02;
            testCase.k = 100;
        end
    end
    
    methods (Test)
        
        function test_a0_coeff(testCase)
            scheme = Ord2ndHelmholtz2D(testCase.k, testCase.h);
            t = {'not expected result for coefficient a0. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.a0, 0, strjoin(t));
        end
                
        function test_as_coeff(testCase)
            scheme = Ord2ndHelmholtz2D(testCase.k, testCase.h);
            t = {'not expected result for coefficient as. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.as, 1, strjoin(t));
        end
                
        function test_ac_coeff(testCase)
            scheme = Ord2ndHelmholtz2D(testCase.k, testCase.h);
            t = {'not expected result for coefficient ac. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.ac, 0, strjoin(t));
        end
                        
    end
    
end

