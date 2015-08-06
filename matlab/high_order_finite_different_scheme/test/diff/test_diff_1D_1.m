% test to show how to differentiate a simple 1D function two time

close all; clear variables; clc;

% possible error formula
mod = @(anal, comp, y) norm(abs(anal - comp), y);
err_rel = @(anal, comp) mod(anal, comp,2) ./ mod(anal, 0, 2); 
norm_err = @(anal, comp) mod(anal, comp,2);
norm_inf = @(anal, comp) mod(anal, comp,inf);
disp_err = @(type_diff, anal, comp) sprintf(strcat(type_diff, ...
    ': \nerr relative (%.2e)', '\nnorm err (%.2e)', '\nnorm inf (%.2e)'), ...
    err_rel(anal, comp), norm_err(anal, comp), norm_inf(anal, comp));

% -------------------------------------------------------------------------
% we define a domain, nb_pt
% -------------------------------------------------------------------------
d = [0,2*pi];
nb_pt = 1000;
h = abs( d(2) - d(1) )/nb_pt ;
x = linspace(d(1), d(2), nb_pt);

% -------------------------------------------------------------------------
% we define a function
% -------------------------------------------------------------------------
func = @(x) sin(x);
func_deriv1 = @(x) cos(x);
func_deriv2 = @(x) -sin(x);
% func = @(x) x.^2 + x + 1;
% func_deriv1 = @(x) 2 * x + 1;

% -------------------------------------------------------------------------
% computation of the elementary differences
% -------------------------------------------------------------------------
deriv_dom = linspace( d(1) -  2*h, d(2) + 2*h, nb_pt + 4 );
fdd = func(deriv_dom);
fx = fdd(3:end-2);
fxph = fdd(4:end-1);
fxmh = fdd(2:end-3);
fxm2h = fdd(1:end-4);

% -------------------------------------------------------------------------
% computation of the first derivative
% -------------------------------------------------------------------------
% forward difference (O(h))
d1f = (fxph - fx)  ./ h;
% backward difference (O(h))
d1b = (fx - fxmh)  ./ h;
% central difference (O(h^2))
d1c = (fxph - fxmh)  ./ (2*h);
d1c2 = (d1f + d1b)  ./ 2;
% order (O(h^3))
d1o3 = (1/(6*h)) * (2 * fxph + 3 * fx - 6 * fxmh + fxm2h);
% analytic derivative
d1a = func_deriv1(x);

disp_err('forward difference', d1a, d1f)
disp_err('backward difference', d1a, d1b)
disp_err('central difference', d1a, d1c)
disp_err('central difference (formula 2)', d1a, d1c2)
disp_err('Order 3 difference', d1a, d1o3)

% -------------------------------------------------------------------------
% computation of the second derivative
% -------------------------------------------------------------------------
% computation of the second derivative
d2c = (fxph - 2 * fx + fxmh) ./ (h.^2);
% analytic second derivative
d2a = func_deriv2(x);
disp_err('second derivative', d2a, d2c)



















