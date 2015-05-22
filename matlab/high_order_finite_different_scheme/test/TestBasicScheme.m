classdef TestBasicScheme < matlab.unittest.TestCase                                
    %TESTBASICSCHEME test the class ORD2NDHEMHOLTZ2D    
    
    methods (Test)
        
        function test_number_line(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 0;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            t = {'not expected result for number of line '...
                'Please recheck the class BasicScheme'};
            testCase.verifyEqual(bs.m, 5, strjoin(t));
        end
        
        function test_number_column(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 0;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            t = {'not expected result for number of column '...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(bs.n, 5, strjoin(t));
        end
                
        function test_central_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 0;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 3; j = 3;
            [c_A, v_A, c_b, v_b] = bs.c_pt(i,j);
            mn = param.m * param.n;
            A = sparse(mn, mn);
            b = sparse(mn);
            A(c_A) = v_A;
            b(c_b) = v_b;
            %  full(A) % debug
            label = j + ( param.m - i ) * param.m;  
            exp_A = [0 0 0 0 0 0 0 1 0 0 0 1 -3 1 0 0 0 1 0 0 0 0 0 0 0];
            exp_b = 0;
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));
            
        end
                        
        function test_north_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 5; j = 4;
            [c_A, v_A, c_b, v_b] = bs.n_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [0 0 1 -3 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            exp_b = -5;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end        

                        
        function test_east_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 2; j = 5;
            [c_A, v_A, c_b, v_b] = bs.e_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 -3 0 0 0 0 1];
            exp_b = -5;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end                
                        
        function test_south_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 1; j = 3;
            [c_A, v_A, c_b, v_b] = bs.s_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 -3 1 0];
            exp_b = -5;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end                     

        function test_west_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 2; j = 1;
            [c_A, v_A, c_b, v_b] = bs.w_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -3 1 0 0 0 1 0 0 0 0];
            exp_b = -5;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end

        function test_north_east_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 5; j = 5;
            [c_A, v_A, c_b, v_b] = bs.ne_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [0 0 0 1 -3 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            exp_b = -10;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end

        function test_south_east_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 1; j = 5;
            [c_A, v_A, c_b, v_b] = bs.se_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 -3];
            exp_b = -10;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end
        
        function test_south_west_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 1; j = 1;
            [c_A, v_A, c_b, v_b] = bs.sw_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -3 1 0 0 0];
            exp_b = -10;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end

        function test_north_west_pt(testCase)
            param.m = 5;
            param.n = 5;            
            param.dirichlet = @(i,j) 5;
            scheme = Ord2ndHelmholtz2D(100, 0.01);
            bs = BasicScheme(param, scheme);
            
            i = 5; j = 1;
            [c_A, v_A, c_b, v_b] = bs.nw_pt(i,j);                        
%             c_A 
%             v_A 
%             c_b 
%             v_b            
            mn = param.m * param.n;            
            A = sparse(mn, mn);
            b = sparse(mn);
            
            A(c_A) = v_A;
            b(c_b) = v_b;
            
%             full(A) % debug
            label = j + ( param.m - i ) * param.m;              
            exp_A = [-3 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            exp_b = -10;
            
            t = {'matrix A, not expected coefficient at line: ', ...
                sprintf('%d',label),...
                'Please recheck the class BasicScheme'};            
            testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
            
            t = {'vector b, not expected coefficient at line: ', ...
                sprintf('%d',label), ...
                'Please recheck the class BasicScheme'};                        
            testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
        end        
        
%         function test_north_pt_sommerfeld(testCase)
%             param.m = 5;
%             param.n = 5;            
%             param.dirichlet = @(i,j) 5;
%             param.k = 100; 
%             param.h = 0.01;
%             
%             param.north = 'sommerfeld';
%             scheme = Ord2ndHelmholtz2D(param.k, param.h);
%             sommerfeld = Ord6thSommerfeld2D( param.h, param.k);
%             bss = BasicSommScheme(scheme, sommerfeld);
%             bs = BasicScheme(param, scheme, bss);            
%             
%             i = 5; j = 4;
%             [c_A, v_A, c_b, v_b] = bs.n_pt(i,j);                        
%        
%             mn = param.m * param.n;            
%             A = sparse(mn, mn);
%             b = sparse(mn);
%             
%             A(c_A) = v_A;
%             b(c_b) = v_b;
%             
%             label = j + ( param.m - i ) * param.m;              
%             exp_A = [0 0 1 (-3.000000000000000 -84.166666666666671i) ...
%                 1 0 0 0 51 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%             exp_b = 0;
%             
%             t = {'matrix A, not expected coefficient at line: ', ...
%                 sprintf('%d',label),...
%                 'Please recheck the class BasicScheme'};            
%             testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
%             
%             t = {'vector b, not expected coefficient at line: ', ...
%                 sprintf('%d',label), ...
%                 'Please recheck the class BasicScheme'};                        
%             testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
%         end        

%         function test_east_pt_sommerfeld(testCase)
%             param.m = 5;
%             param.n = 5;            
%             param.dirichlet = @(i,j) 5;
%             param.k = 100; 
%             param.h = 0.01;
%             
%             param.east = 'sommerfeld';
%             scheme = Ord2ndHelmholtz2D(param.k, param.h);
%             sommerfeld = Ord6thSommerfeld2D( param.h, param.k);
%             bss = BasicSommScheme(scheme, sommerfeld);
%             bs = BasicScheme(param, scheme, bss);            
%             
%             i = 4; j = 5;
%             [c_A, v_A, c_b, v_b] = bs.e_pt(i,j);                        
%        
%             mn = param.m * param.n;            
%             A = sparse(mn, mn);
%             b = sparse(mn);
%             
%             A(c_A) = v_A;
%             b(c_b) = v_b;
%             
%             label = j + ( param.m - i ) * param.m;              
%             exp_A = [0 0 1 (-3.000000000000000 -84.166666666666671i) ...
%                 1 0 0 0 51 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%             exp_b = 0;
%             
%             t = {'matrix A, not expected coefficient at line: ', ...
%                 sprintf('%d',label),...
%                 'Please recheck the class BasicScheme'};            
%             testCase.verifyEqual(full(A(label,:)), exp_A, strjoin(t));
%             
%             t = {'vector b, not expected coefficient at line: ', ...
%                 sprintf('%d',label), ...
%                 'Please recheck the class BasicScheme'};                        
%             testCase.verifyEqual(full(b(label)), exp_b, strjoin(t));            
%         end         
        
    end
    
end

