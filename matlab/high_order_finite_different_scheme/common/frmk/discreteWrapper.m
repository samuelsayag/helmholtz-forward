function [ wrapper ] = discreteWrapper( x1, x2, h, wrappee )
%INDXFUNCWRAPPER Being given a 2D domain [x1,y1]*[x2,y2] (and a step h) on
%which the function wrappee has a range, we build a wrapper that gives the
%value of the function receiving a discrete index.
%wrappee(x1,x2) = wrapper(1,1) and wrappee(x1+h,x2) = wrapper(2,1) and so
%on...

x_val = @(x) x1 + (x-1) * h;
y_val = @(x) x2 + (x-1) * h;

wrapper = @(i,j) wrappee( x_val(i), y_val(j) );

end

