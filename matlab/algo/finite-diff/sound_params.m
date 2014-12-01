function [ c, omega, k, k2 ] = sound_params()
% function [ s, f, k, k2 ] = sound_params()
% Produce sound parameters for finite difference
% c : vector of speed of sound in tissue 
% f : frequency
% k : k = (2 * pi * f) / c


% s = 1450:25:1650
c = 1400:50:1700;

f = 0.5e6:0.5e6:20e6;

omega = 2 * pi * f ; 

k = omega' * (1./ c);

k2 = power(k,2);

end

