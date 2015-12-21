function [ nb, h ] = GridSize( k, so )
%GRIDSIZE provide an evaluation of the number of point of a grid with
%parameter:
% k: the wave number
% so: scheme order (scheme polynom error order)

nb = k * ( (so+1) / so );
h = 1/nb;

end

