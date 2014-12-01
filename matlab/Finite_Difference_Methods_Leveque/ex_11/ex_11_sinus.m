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
% the "real" function
du_v = cos(x0);

format long

% display error
% Note: the order of the error
back_err = du_v - du_moins; 
forw_err = du_v - du_plus;
cent_err = du_v - du_cent;
estim_err = [h, back_err, forw_err, cent_err];
estim_err


% graph the err = f(h) 
plot(log(h), log(abs(back_err)), '.-', log(h), log(forw_err), '.-', log(h), log(cent_err), '.-')
% plot(h, back_err, h, forw_err, h, cent_err, '.-')
title('Error (back, forw, cent) against step (h)')
xlabel('step (h)')
ylabel('errors')
legend('backward', 'forward', 'central', 'Location', 'NorthWest')

format longg