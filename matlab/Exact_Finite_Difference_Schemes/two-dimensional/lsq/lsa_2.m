function [ res ] = lsa_1( p, af, lsp )
%LSA_1 - Function that compute the optimal angle theta necessary for a
%simulation of the Helmholtz equation being given:
% param:
%   - p: structure of parameter necessary for a simulation of the Helmholtz
%   equation.
%   - af: the function that give the analytic solution of the Helmholtz
%   equation.
%   - lsp: the parameter of the least square algorithm. An structure of 4
%   components [ theta1, theta2, na, pr ]. The two first parameters theta1 
%   and theta2 give the range of angle between them is computed the Bessel
%   function for the central point coefficient. The parameter na is the
%   number of angle inside this range that let us yields the different
%   angles to initialize the non linear least-square search of the optimal
%   theta as well as compute the sommerfeld constraint. 
%   The parameter pr is the precision that we want to reach on theta. 

check_parameters(af, lsp);
res = {};
% t1 and t2 are the angle that define the range we wish to narrow
t1 = min([lsp.theta1 lsp.theta2]);  
t2 = max([lsp.theta1 lsp.theta2]);  
% these are the bounds we used to drop the optimized angle that are not
% satisfying us from a range point of view
t_min = t1; t_max = t2;
% buid a factice condition just to enter the loop
converge = lsp.pr + 1;

while converge > lsp.pr  
    
    % STEP 2. - Compute the angles between theta1 and theta2 => ra
    ra = (linspace( 0, lsp.na, lsp.na+1 ) * (t2 - t1)/lsp.na + t1)';
    c_ra = mat2cell( ra , ones(size(ra)));
    
    % STEP 2. - Compute the angles between theta1 and theta2 => ra
    bsp = @(t) build_sim_param(t, t1, t2, p);
    ps = cellfun(bsp, c_ra, 'UniformOutput', false);
    
    sh = @(p) solve(p);
    [As, bs, xs] =  cellfun(sh, ps, 'UniformOutput', false);
    
    % STEP 2. - Compute F cell of all f such that: f = x(j) - x(theta,j)
    fn = @(x) @(th) norm(x - af(th)) ; % the error function to minimize ( = f )    
    Fns = cellfun(fn, xs, 'UniformOutput', false);
    
    % STEP 2. - Compute the non linear least square solution for F (all f's)
    ln = @(fn, x0) lsqnonlin(fn, x0); % lsq algo we wish to apply for all theta    
    [x_n, resnorm_n, residual_n, exitflag_n ] = cellfun(ln, Fns, c_ra, 'UniformOutput', false);

    x = x_n;    
    res = [res ; cell2mat([c_ra x, resnorm_n, exitflag_n])];
    
    % STEP 2.4 - Extract the angles for them solution is between the set
    % theta1 and theta2. Get the new theta1 and theta2 (min, max)(ag_ext)
    x = cell2mat(x);
    in_range = x(and(x < t_max,  x > t_min));
    t1 = min(in_range); t2 = max(in_range);
    [t1, t2]
    % STEP 2.3 - Compute the distance between theta1 and theta2. If the
    % distance is satisfying the precision exit.
    converge = t2 - t1;
end
   
end 

function [ p ]  = build_sim_param(theta, t1, t2, p)
    p.theta = theta;
    p.bessel = @(x) bessel_integral(x, t1, t2);    
end

function [ A, b, x ] = solve( p )
    [ func_scheme, p ] = helmholtz_2D_scheme_factory( p );
    [ A, b ] = build_two_dimensional_problem2(p, func_scheme);
    sol = A\b;    
    % tmp_sol = bicgstab(A, b, 1e-6, 1e5);
    % tmp_sol = qmr(A, b, 1e-6, 1e5);
    % tmp_sol = gmres(A, b, 60, 1e-6, 300);
    % important to leave the transpose function for it is a complex matrice
    x = transpose(reshape(sol, p.m, p.n));
end

function [af, lsp] = check_parameters(af, lsp)
    display(sprintf('check the parameters...'));
end