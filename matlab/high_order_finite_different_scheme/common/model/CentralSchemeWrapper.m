classdef CentralSchemeWrapper < NinePtStencil
    %SCHEMEWRAPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        scheme = [];
    end
    
    methods (Access = public)
        
        function obj = CentralSchemeWrapper(scheme)
            % this function receive a CentralScheme interface and check it
            obj.check_param(scheme);
            obj.scheme = scheme;
        end         
        
        % the central point coefficient
        function c = c(obj)
            c = obj.scheme.a0;
        end
        % the north point coefficient
        function n = n(obj)
            n = obj.scheme.as;
        end
        % the north east point coefficient
        function ne = ne(obj)
            ne = obj.scheme.ac;
        end
        % the east point point coefficient
        function e = e(obj)
            e = obj.scheme.as;
        end
        % the south east point coefficient
        function se = se(obj)
            se = obj.scheme.ac;
        end
        % the south point coefficient
        function s = s(obj)
            s = obj.scheme.as;
        end
        % the south west point coefficient
        function sw = sw(obj)
            sw = obj.scheme.ac;
        end
        % the west point coefficient
        function w = w(obj)
            w = obj.scheme.as;
        end
        % the north west point coefficient
        function nw = nw(obj)      
            nw = obj.scheme.ac;
        end
    end
    
    methods (Access = private)
        function obj = check_param(obj, scheme)
            narginchk(2,2);
            p = inputParser;
            addRequired(p, 'scheme', @(s) isa(s,'CentralScheme'));
            parse(p, scheme);
        end
    end    
    
end

