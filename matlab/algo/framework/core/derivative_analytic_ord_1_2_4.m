function [ f, fx, fy, fxx, fyy, fxxyy, fxxxx, fyyyy ] = ...
    derivative_analytic_ord_1_2_4( clos_f )
%DERIVATIVE_ANALYTIC compute a closure (pointer to function) of the
%function f and its first and second third and fourth derivative 
%with respect to x and y
%
%clos_f: a closure of the mathematical expression of the function ex: 
%x.^2+x+1 or sin(x) + log(y)
%f, fx, fy, fxx, fyy, fxxyy, fxxxx, fyyyy : closure of the function f : 
% (@(x,y) f(x,y)) and closure to all of the derivative.
%
%Usage: 
%> f = @(x,y) sin(x + y );
%> [ f, fx, fy, fxx, fyy, fxxyy, fxxxx, fyyyy ] = derivative_analytic( f );
%

syms x y; % declare 2 symbolic variables
tof = @( z ) matlabFunction( z , 'Vars', [x,y]);

fs = sym(clos_f);
f = tof(fs);

fxs = diff( fs, x );
fx = tof(fxs);

fys = diff( fs, y );
fy = tof(fys);

fxxs = diff(fxs, x);
fxx = tof(fxxs);

fyys = diff(fys, y);
fyy = tof(fyys);

fxxys = diff(fxxs, y);
fxxyys = diff(fxxys, y);
fxxyy = tof(fxxyys);

fxxxs = diff(fxxs, x);
fxxxxs = diff(fxxxs, x);
fxxxx = tof(fxxxxs);

fyyys = diff(fyys, y);
fyyyys = diff(fyyys, y);
fyyyy = tof(fyyyys);

end

