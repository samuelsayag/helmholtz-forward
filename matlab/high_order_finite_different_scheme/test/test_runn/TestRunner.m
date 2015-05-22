clear all; clc;

import matlab.unittest.TestSuite;

a = TestSuite.fromFile('TestBasicScheme.m');
b = TestSuite.fromFile('TestDirichlet.m');
c = TestSuite.fromFile('TestExactHemholtz2D.m');
d = TestSuite.fromFile('TestExactSommerfeld2D.m');
e = TestSuite.fromFile('TestMatrixBuilder.m');
f = TestSuite.fromFile('TestOrd2ndHemholtz2D.m');
g = TestSuite.fromFile('TestOrd2ndSommerfeld2D.m');
h = TestSuite.fromFile('TestOrd4thHemholtz2D.m');
i = TestSuite.fromFile('TestOrd6thHemholtz2D.m');
j = TestSuite.fromFile('TestOrd6thSommerfeld2D.m');

largeSuite = [a, b, c, d, e, f, g, h, i, j];

result = run(largeSuite);