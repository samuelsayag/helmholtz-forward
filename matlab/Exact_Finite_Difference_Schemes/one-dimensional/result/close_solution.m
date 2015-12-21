clear all;

% define the close probleme domain and parameters
a0 = 0; 
a1 = 1;
nbr_pt = 101;
x = linspace(a0, a1, nbr_pt);
k = [10 20 50 100 150];

figure(1); 
cpt_fig = 1;
for i = k
    % get the close solution
    y = analytic_sol_1D(i, x);
    
    % display the real part
    subplot(size(k,2),1, cpt_fig);
    plot(x, real(y));
    title(sprintf('Real part, k = %u, ', i));
    cpt_fig = cpt_fig + 1;
    
%     % display the imaginary part
%     subplot(size(k,2),2, cpt_fig);
%     plot(x, imag(y));
%     title(sprintf('Imaginary part, k = %u, ', i));
%     cpt_fig = cpt_fig + 1;    
end

