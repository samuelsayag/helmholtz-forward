classdef TestOrd6thHemholtz2D < matlab.unittest.TestCase                                
    %TESTORD6THHEMHOLTZ2D test the class ORD2NDHEMHOLTZ2D    
    properties
        h;
        k;
    end
 
%     methods(TestMethodSetup)
%         function create_h_k(testCase)
%             testCase.h = 0.01;
%             testCase.k = 100;
%         end
%     end
    
    methods (Test)
        
        function test_a0_coeff_1(testCase)
            scheme = Ord6thHelmholtz2D(sqrt(30), 1, 3);
            t = {'not expected result for coefficient a0. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.a0, 57/3, strjoin(t));
        end
        
        function test_a0_coeff_2(testCase)
            scheme = Ord6thHelmholtz2D(100, 0.01, 5);
            testCase.verifyEqual(scheme.delta, 5, 'gamma not correct');
            t = {'not expected result for coefficient a0. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.a0, 68/90 - 10/3, strjoin(t));
        end
        
        function test_as_coeff_1(testCase)
            scheme = Ord6thHelmholtz2D(100, 0.01, 1.5);
            t = {'not expected result for coefficient as. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.as, 2/3 + 2/45, strjoin(t));
        end        
        
        function test_as_coeff_2(testCase)
            scheme = Ord6thHelmholtz2D(100, 0.01, -13/2);
            t = {'not expected result for coefficient as. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.as, 2/3 + 3/45, strjoin(t));
        end        
         
        function test_ac_coeff_1(testCase)
            scheme = Ord6thHelmholtz2D(100, 0.01, 0);
            t = {'not expected result for coefficient ac. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.ac, 1/6 + 7/360, strjoin(t));
        end
         
        function test_ac_coeff_2(testCase)
            scheme = Ord6thHelmholtz2D(100, 0.01, 2);
            t = {'not expected result for coefficient ac. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.ac, 1/6 + 8/360, strjoin(t));
        end
                        
        function test_bs(testCase)
            scheme = Ord6thHelmholtz2D(testCase.k, testCase.h);
            t = {'not expected result for coefficient bs. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.bs, 1, strjoin(t));
        end        
                        
        function test_bc(testCase)
            scheme = Ord6thHelmholtz2D(testCase.k, testCase.h);
            t = {'not expected result for coefficient bc. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.bc, 0, strjoin(t));
        end                        
    end
    
end

