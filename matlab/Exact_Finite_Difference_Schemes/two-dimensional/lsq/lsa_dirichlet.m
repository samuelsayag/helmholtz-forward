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
t1 = min([lsp.theta1 lsp.theta2]); 
t_min = t1;
t2 = max([lsp.theta1 lsp.theta2]); 
t_max = t2;
converge = lsp.pr + 1;

while converge > lsp.pr 
    % STEP 1 - Build the system Ax=b and solve it => x_t (temporary solution)
    p.bessel = @(x) bessel_integral(x, t1, t2);
    [ A, b, x_t ] = solve( p );    
    
    % STEP 2.1 - Compute the angles between theta1 and theta2 => ra
    ra = linspace( 0, lsp.na, lsp.na+1 ) * (t_max - t_min)/lsp.na + t_min;
%     ra = linspace( 0, lsp.na, lsp.na+1 ) * (t2 - t1)/lsp.na + t1;
    c_ra = mat2cell( ra' , ones(size(ra)));

    % STEP 2.2 - Compute F matrix of all f such that: f = x_t(j) - x(theta,j)
    sample_f = @(a) [a(:,1); a(end, 2:end)'; a(1:size(a,1)-1, 2); a(size(a,1)-1, 3:end)'];
%     sample_f = @(a) [a(:,1); a(end, 2:end)'];
    fn = @(th) norm(sample_f(x_t) - sample_f(af(th))) ; % the error function to minimize ( = f )
    ln = @(x0) lsqnonlin(fn, x0); % lsq algo we wish to apply for all theta    

    % STEP 2.3 - Compute the non linear least square solution for F (all f's)
    [x_n, resnorm_n, residual_n, exitflag_n ] = cellfun(ln, c_ra, 'UniformOutput', false);

    x = x_n;    
    res = [res ; cell2mat([c_ra x, resnorm_n, exitflag_n])];
    
    % STEP 2.4 - Extract the angles for them solution is between the set
    % theta1 and theta2. Get the new theta1 and theta2 (min, max)(ag_ext)
    x = cell2mat(x);
    x'
    in_range = x(and(x <= t2,  x >= t1));
    t1 = min(in_range); t2 = max(in_range);
    [t1, t2]
    % STEP 2.3 - Compute the distance between theta1 and theta2. If the
    % distance is satisfying the precision exit.
    converge = t2 - t1;
end
   
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