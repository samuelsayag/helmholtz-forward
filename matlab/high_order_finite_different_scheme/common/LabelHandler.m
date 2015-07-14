function [ label ] = LabelHandler( nbline, nbcol, i, j )
%LABELHANDLER compute the standard label for a line of the A and b vector
% matrix to solve in the context of a finite difference system ( A x = b )
% so that the resulting matrix is k-diagonal.

    label = j + (nbline - i) * nbcol;

end

