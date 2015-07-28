classdef (Abstract) NinePtStencil
    %NINEPTSTENCIL this abstract definition state the methods that a nine
    % point scheme must expose to be processed by the BasicScheme.
    % 
    % This stencil assume coefficient for this 
    %
    %     nw ------ n ------ ne
    %      |        |        |
    %      w ------ c ------ e  
    %      |        |        |
    %     sw ------ s ------ se
    %
    methods (Abstract)
        % the central point coefficient
        c = c(obj);
        % the north point coefficient
        n = n(obj);
        % the north east point coefficient
        ne = ne(obj);
        % the east point point coefficient
        e = e(obj);
        % the south east point coefficient
        se = se(obj);
        % the south point coefficient
        s = s(obj);
        % the south west point coefficient
        sw = sw(obj);
        % the west point coefficient
        w = w(obj);
        % the north west point coefficient
        nw = nw(obj);
    end
    
end

