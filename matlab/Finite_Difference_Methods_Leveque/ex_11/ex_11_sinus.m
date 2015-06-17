close all; clear variables; clc;
% the step size we want to try
h = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3]';

x0 = ones(5,1);
u = sin(x0);

% the forward difference
du_plus = (sin(x0 + h) - sin(x0))./h;
% the backward difference
du_moins = (sin(x0) - sin(x0-h))./h;
% the central difference
du_cent = 1/2 * (du_plus + du_moins);
% the O(h^3)
d3_u = 1./(6 * h) .* (2 * sin(x0 + h) + 3 * sin(x0) - 6 * sin(x0 - h) + sin(x0 - 2 * h));  
% the "real" function
du_v = cos(x0);

format long

% display error
% Note: the order of the error
back_err = du_v - du_moins; 
forw_err = du_v - du_plus;
cent_err = du_v - du_cent;
du_3_err = du_v - d3_u;
res = {'h' 'Back Err' 'Forw Err' 'Cent Err' 'du_3_err'};
f = @(m) mat2cell(m,ones(5,1));
res = [res; f(h) f(back_err) f(forw_err) f(cent_err) f(du_3_err)];
res


% graph the err = f(h) 
plot(log(h), log(abs(back_err)), '.-', log(h), log(forw_err), '.-', ...
    log(h), log(cent_err), '.-', log(h), log(abs(du_3_err)), '.-')
% plot(h, back_err, h, forw_err, h, cent_err, '.-')
title('Error (back, forw, cent) against step (h), log|E.h| = log(C) + p * log(h)')
xlabel('step (h)')
ylabel('errors')
legend('backward', 'forward', 'central', 'du_3 O(h^3)', 'Location', 'SouthEast')

format longg