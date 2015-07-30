function [ wrapper ] = indxFuncWrapper( x1, x2, h, wrappee )
%INDXFUNCWRAPPER wrapp the handle to and index dependent one

x_val = @(x) x1 + (x-1) * h;
y_val = @(x) x2 + (x-1) * h;

wrapper = @(i,j) wrappee( x_val(i), y_val(j) );


end

