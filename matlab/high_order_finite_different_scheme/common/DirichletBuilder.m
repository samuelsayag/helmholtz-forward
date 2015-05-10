function [ discrete_dirichlet ] = DirichletBuilder( param )
%DIRICHLETBUILDER This function build a function handle for the dirichlet 
% boundary. It takes as a parameter a function of the type f(x,y) with x
% and y being the coordinate and transforms it to a function f(i,j) where i
% and j are the indices of each point of the grid on which is simulated the
% problem with a step h.
% the parameter param is assumed to be a struct that contain at least the
% following fields a, c, h. The area of the grid being simulated for x
% belonging to [ a, b ] and y belong to [ c, d ] and h is the step between
% to successive indices i and i+1.

x_val = @(x) param.a + (x-1) * param.h;
y_val = @(x) param.c + (x-1) * param.h;

discrete_dirichlet = ...
    @(i,j) param.dirichlet( x_val(j), y_val(i) );

end

