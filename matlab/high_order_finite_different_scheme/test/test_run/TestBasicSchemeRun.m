close all; clearvars; clc;

suite = TestBasicScheme();
res = run(suite);
res

% import matlab.unittest.TestSuite;
% suite   = TestSuite.fromMethod(?TestBasicScheme, ...
%     'test_north_pt_sommerfeld');


% import matlab.unittest.TestSuite;
% suite   = TestSuite.fromMethod(?TestBasicScheme, 'test_central_pt');
% res = run(suite);



