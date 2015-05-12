classdef TestOrd4thHemholtz2D < matlab.unittest.TestCase                                
    %TESTORD2NDHEMHOLTZ2D test the class ORD2NDHEMHOLTZ2D    
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
            scheme = Ord4thHelmholtz2D(sqrt(2) * 100, 0.01);
            t = {'not expected result for coefficient a0. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.a0, -2, strjoin(t));
        end
        
        function test_a0_coeff_2(testCase)
            scheme = Ord4thHelmholtz2D(sqrt(2) * 100, 0.01, 18);
            testCase.verifyEqual(scheme.gamma, 18, 'gamma not correct');
            t = {'not expected result for coefficient a0. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.a0, -1, strjoin(t));
        end
        
        function test_as_coeff_1(testCase)
            scheme = Ord4thHelmholtz2D(200, 0.01);
            t = {'not expected result for coefficient as. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.as, 1, strjoin(t));
        end        
        
        function test_as_coeff_2(testCase)
            scheme = Ord4thHelmholtz2D(200, 0.01, 72);
            t = {'not expected result for coefficient as. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.as, -3, strjoin(t));
        end        
         
        function test_ac_coeff(testCase)
            scheme = Ord4thHelmholtz2D(100, 0.01, 144);
            t = {'not expected result for coefficient ac. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.ac, 1 + 1/6, strjoin(t));
        end
                        
        function test_bs(testCase)
            scheme = Ord4thHelmholtz2D(testCase.k, testCase.h);
            t = {'not expected result for coefficient bs. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.bs, 1, strjoin(t));
        end        
                        
        function test_bc(testCase)
            scheme = Ord4thHelmholtz2D(testCase.k, testCase.h);
            t = {'not expected result for coefficient bc. '...
                'Please recheck the class Ord2ndHelmholtz2D'};
            testCase.verifyEqual(scheme.bc, 0, strjoin(t));
        end                        
    end
    
end

