function [A, b, x, params] = intStd_diriW_diriE()
% std_diriW_diriE
% Summary of this function goes here
%   Detailed explanation goes here

close all; clear all; clc;

x = [];

% parameters necessary to compute interior points
params.k = 10;
params.h = 0.5./params.k;
a = 0;
b = 1;
params.m = (b-a)./params.h;
params.interior = 'std';
params.boundary = 'dirichlet';
params.dirichlet.W = @(params, A, b, i) 1;
params.dirichlet.E = @(params, A, b, i) 1;

% params

[ func_scheme, params ] = helmholtz_1D_scheme_factory( params );

% create the matrix of finite difference
[A, b] = build_one_dimensional_problem2(params, func_scheme);

% debug
% full(A)
% b

%------------------- solve the system -----------------------------
tstart = tic;

tol = 1e-6;
[x,flag,relres] = gmres(A, b', [], tol);        

telapsed = toc(tstart);
%------------------- display some result -----------------------------

telapsed
