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
dx = [0,2*pi];
dy = [0,2*pi];
nb_pt = 100;
h = abs( dx(2) - dx(1) )/nb_pt ;
x = linspace(dx(1), dx(2), nb_pt);
y = linspace(dy(2), dy(1), nb_pt);
[X, Y] = meshgrid(x,y);

% -------------------------------------------------------------------------
% we define a function and its derivative (to compare with computed
% results)
% -------------------------------------------------------------------------
strf = 'sin( x+y )';
syms x y; % declare 2 symbolic variables
tof = @( z ) matlabFunction( z , 'Vars', [x,y]);

sf = sym(strf);
f_a = tof(sf);

dfxs = diff( sf, x );
dfx_a = tof(dfxs);

dfys = diff( sf, y );
dfy_a = tof(dfys);

dfxxs = diff(dfxs);
dfxx_a = tof(dfxxs);

dfyys = diff(dfys);
dfyy_a = tof(dfyys);

% -------------------------------------------------------------------------
% computation of the elementary differences
% -------------------------------------------------------------------------
dx = linspace(dx(1)-h, dx(2)+h, nb_pt+2);
dy = linspace(dy(2)+h, dy(1)-h, nb_pt+2);
[dX, dY] = meshgrid(dx,dy);

fdd = f_a(dX, dY);

fx = fdd(2:end-1, 2:end-1);
fxph = fdd(2:end-1, 3:end);
fxmh = fdd(2:end-1, 1:end-2);

fy = fdd(2:end-1, 2:end-1);
fyph = fdd(1:end-2, 2:end-1);
fymh = fdd(3:end, 2:end-1);

% -------------------------------------------------------------------------
% computation of the first derivative with respect to X
% -------------------------------------------------------------------------
% forward difference (O(h))
d1xf = (fxph - fx)  ./ h;
% backward difference (O(h))
d1xb = (fx - fxmh)  ./ h;
% central difference (O(h^2))
d1xc = (fxph - fxmh)  ./ (2*h);
d1xc2 = (d1xf + d1xb)  ./ 2;
% analytic derivative
d1xa = dfx_a(X,Y);

disp_err('forward difference', d1xa, d1xf)
disp_err('backward difference', d1xa, d1xb)
disp_err('central difference', d1xa, d1xc)
disp_err('central difference (formula 2)', d1xa, d1xc2)

% -------------------------------------------------------------------------
% computation of the first derivative with respect to X
% -------------------------------------------------------------------------
% forward difference (O(h))
d1yf = (fyph - fy)  ./ h;
% backward difference (O(h))
d1yb = (fy - fymh)  ./ h;
% central difference (O(h^2))
d1yc = (fyph - fymh)  ./ (2*h);
d1yc2 = (d1yf + d1yb)  ./ 2;
% analytic derivative
d1ya = dfy_a(X,Y);

disp_err('forward difference', d1ya, d1yf)
disp_err('backward difference', d1ya, d1yb)
disp_err('central difference', d1ya, d1yc)
disp_err('central difference (formula 2)', d1ya, d1yc2)

% -------------------------------------------------------------------------
% computation of the second derivative with respect to X
% -------------------------------------------------------------------------
% computation of the second derivative
d2xc = (fxph - 2 * fx + fxmh) ./ (h.^2);
% analytic second derivative
d2xa = dfxx_a(X, Y);
disp_err('second derivative (X)', d2xa, d2xc)

% -------------------------------------------------------------------------
% computation of the second derivative with respect to Y
% -------------------------------------------------------------------------
% computation of the second derivative
d2yc = (fyph - 2 * fy + fymh) ./ (h.^2);
% analytic second derivative
d2ya = dfyy_a(X, Y);
disp_err('second derivative (Y)', d2ya, d2yc)




















