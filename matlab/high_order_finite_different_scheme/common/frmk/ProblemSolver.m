classdef ProblemSolver
    %PROBLEMSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        % A scheme of the type:
        % {'Ord2ndHelmholtz2D', 'Ord4thHelmholtz2D', 'Ord6thHelmholtz2D'}
        scheme = [];
        % A scheme for the sommerfeld boundary {'Ord6thSommerfeld2D'}
        sommerfeld = [];
        % A struct that fields are as follow:
        % param.h: length of the step in the unit chosen ex: 0.002 = 2mm if
        %   the meter is chosen as a unit.
        % param.k: the wave number k = 2 * pi * f / c (f = frequency,
        %   c = celerity of the wave)
        % param.m: the number of point along the y axis (number of line)
        % param.n: the number of point along the x axis (number of column)
        % param.dirichlet: the dirichlet function as pointer of function
        %   that give a value by passing it the line and column index. Ex:
        %   dirichlet = @(i,j) some_function....
        % param.north = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the north side of the area.
        % param.east = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the east side of the area.
        % param.south = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the south side of the area.
        % param.west = {'dirichlet', 'sommerfeld'} type of the boundary on
        %   the west side of the area.
        % 
        % Note: it will not be possible to choose Sommerfeld type
        % boundaries along all the side. At least one must be of Dirichlet
        % type.
        param = [];        
        % A solver for the equation A x = b. With A being a matrice and b a
        % vector. It is given as a function handle.
        solver = [];
    end
    
    methods (Access = public)
        
        function obj = ProblemSolver(param, scheme, solver, sommerfeld)            
            narginchk(3, 4);
            obj = obj.check_param(param, scheme);            
            obj.scheme = scheme;
            obj.param = param;
            obj.param.dirichlet = DirichletBuilder( param ); 
            obj = obj.set_solver(solver);
            if nargin == 4
                obj = obj.check_sommerfeld(sommerfeld);
                obj.sommerfeld = sommerfeld;
            end
        end
        
        function obj = set_solver(obj, solver)
           obj = obj.check_solver(solver);
           obj.solver = solver; 
        end
        
        function [A, b, x] = solve(obj)
            % function [ A, b, x ] = solve()
            % This function is responsible to solve once the problem
            % proposed under the form of a param struct that is described
            % in the properties and the given scheme also described in the
            % properties of this class.
            if ~isempty(obj.sommerfeld)
                bss = BasicSommScheme( obj.scheme, obj.sommerfeld );
                bs = BasicScheme(obj.param, obj.scheme, bss);
            else
                bs = BasicScheme(obj.param, obj.scheme);
            end
            mb = MatrixBuilder(bs);
            
            tic
            display(sprintf('build the matrice...'));
            [ A, b ] = mb.build();
            toc
            tic
            display(sprintf('solve Ax=b...'));            
            x = obj.solver(A, b);
            toc
            x = transpose(reshape(x, obj.param.n, obj.param.m));
        end
    end
    
    methods(Access = private)
        function obj = check_param(obj, param, scheme)
            p = inputParser;

            schemes = {'Ord2ndHelmholtz2D',... 
                'Ord4thHelmholtz2D','Ord4thHelmholtz2D_2',...
                'Ord6thHelmholtz2D', 'Ord6thHelmholtz2D_2', ...
                'ExactScheme2D', 'Poisson2D'};
            addRequired(p, 'scheme', ...
                @(x)validateattributes( x, schemes, {'nonempty'}));
                        
            function res = valide_param(x)
                validateattributes( x, {'struct'}, {'nonempty'});
                res = isfield(x, 'dirichlet');
                validateattributes( x.dirichlet, {'function_handle'}, {'nonempty'});
                res = isfield(x, 'm') && res;
                validateattributes( x.m, {'double'}, {'nonempty', 'integer'});                
                res = isfield(x, 'n') && res;
                validateattributes( x.n, {'double'}, {'nonempty','integer'});
            end
            addRequired(p, 'param', @(x)valide_param(x));                 
            
            parse(p, scheme, param);            
        end
      
        function obj = check_sommerfeld(obj, sommerfeld)
            p = inputParser;            
            sommerfelds = {'Ord6thSommerfeld2D', 'Ord2ndSommerfeld2D', ...
                'ExactSommerfeld2D'};
            addRequired(p, 'sommerfeld', ...
                @(x)validateattributes( x, sommerfelds, {'nonempty'}));        
            parse(p, sommerfeld);           
        end
        
        function obj = check_solver(obj, solver)
            p = inputParser;
            
           function valid_handle(x)
                validateattributes( x, {'function_handle'}, {'nonempty'});
            end
            addRequired(p, 'solver', @(x)valid_handle(x))            
            
            parse(p, solver);
        end    
    end
    
end

