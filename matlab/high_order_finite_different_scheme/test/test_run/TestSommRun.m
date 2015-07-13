clear variables; close all; clc;

suite = TestSommerfeld();

% import matlab.unittest.TestSuite;
% suite   = TestSuite.fromMethod(?TestBasicScheme, ...
%     'test_north_pt_sommerfeld');

res = run(suite);
res