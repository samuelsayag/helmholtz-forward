function [ R_2 ] = linear_reg( source, fitted )
%LINEAR_REG 
% source: the data of the model
% fit: the fitted data
%
% from: 
% http://www.mathworks.com/help/releases/R2015a/matlab/data_analysis/linear-regression.html

% Compute the residual values as a vector of signed numbers:
yresid = source - fitted;
% Square the residuals and total them to obtain the residual sum of squares:
SSresid = sum(yresid.^2);
% Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(source)-1) * var(source);
% Compute R2 using the formula given in the introduction of this topic:
R_2 = 1 - SSresid/SStotal;

end

