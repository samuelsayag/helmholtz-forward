function [ func_scheme, params ] = helmholtz_1D_scheme_factory( params )
%HELMHOLTZ_1D_SCHEME_FACTORY 
% The purpose of this function is to offer a flexible way of creating a
% scheme of the Helmholtz equation. The scheme is computed by considering a
% linear region. One side is a source (the computation of the matrix
% will be similar as a Dirichlet Boundary Condition) but the function of
% the source is left to the client. The three other sides may be
% constrained by Dirichlet Boundary Condition or Sommerfeld condition.
%
% return: 
%   func_scheme: a function that may directly be used as a scheme to create 
% the matrix.n
% parameters:
%   param: a struct that contains
%       m: number of point 
%       k: omega/c = 2*pi / c*T = 2*pi*f / c = 2*pi/lambda 
% (lambda = wavelength, f = frequency, T = period, c = celerity of the wave)
%       h: step choosen for the scheme
%       interior: 'new', 'std'
%       boundary: 'sommerfeld_new', 'sommerfeld_std', 'dirichlet'
% (if dirichlet is set all the side are considered of Dirichlet kind)
%       dirichlet.W: (depend of boundary) function handle
%       dirichlet.E: (depend of boundary) function handle
%       Note:
%           - a Dirichlet function is of signature: f(A, b, i, j) with 'A'
% the scheme matrix ((m*n)²), 'b' the vector such that the solution x of 
% A x = b is the solution of the Helmholtz equation.
%           - AT LEAST ONE of the side must be of Dirichlet Boundary kind.
%           - if a Dirichlet function is not given the side is considered 
% under a Sommerfeld Boundary.

% test if the south_source parameter is a well defined function
check_params(params)

dirich_w = isfield(params.dirichlet, 'W'); 
dirich_e = isfield(params.dirichlet, 'E'); 

schemes.interior_func = factory_central_scheme(params);

key = {'W', dirich_w, params.interior, params.boundary};
schemes.west_bound =  factory_side_scheme(key, params);

key = {'E', dirich_e, params.interior, params.boundary};
schemes.east_bound =  factory_side_scheme(key, params);

func_scheme = @(params, A, b, i, j) generic_scheme_creation(...
    params, A, b, i, j, schemes);

end

function [f_ctrl_pt] = factory_side_scheme(key, params)
    % define small function handle that are valids for all the function
    dirich_w = @(params, A, b, i, j) dirichlet_wrapper(params,...
        A, b, i, j, params.dirichlet.W);
    dirich_e = @(params, A, b, i, j) dirichlet_wrapper(params,...
        A, b, i, j, params.dirichlet.E);
    
    c_std = @(params, A, b, i, j)ctrl_pt_std(params, A, b, i, j);
    c_new = @(params, A, b, i, j)ctrl_pt_new(params, A, b, i, j);
    %dirichlet    
    % define an equal operator for dirichlet choice
    condition = key(1:3);
    ceq = @(x,y) strcmp(x{1}, y{1}) && (x{2} == y{2}) && strcmp(x{3},y{3});
    if (ceq(condition, {'S', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j) dirichlet_side_generic(...
            params, A, b, i, j, c_std, dirich_s,...
            'up_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'S', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, dirich_s,...
            'up_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'W', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_std, dirich_w,...
            'down_pt', 'up_pt', 'right_pt');
        
    elseif (ceq(condition, {'W', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, dirich_w,...
            'down_pt', 'up_pt', 'right_pt');
        
    elseif (ceq(condition, {'N', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_std, dirich_n,...
            'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'N', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, dirich_n,...
            'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'E', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_std, dirich_e,...
            'down_pt', 'left_pt', 'up_pt');
        
    elseif (ceq(condition, {'E', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, dirich_e,...
            'down_pt', 'left_pt', 'up_pt');
    end
    
    % sommerfeld       
    % re-define an equal operator for sommerfeld choice
    condition = key([1,2,4]);
    
    if (ceq(condition, {'S', 0, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_std(...
            params, A, b, i, j, 1, 'up_pt', 'right_pt', 'left_pt');
        
    elseif (ceq(condition, {'S', 0, 'new'}))
        k2h = params.k * sin(params.theta);
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, b, i, j, 1, k2h, 'up_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'W', 0, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_std(...
            params, A, b, i, j, 1, 'right_pt', 'down_pt', 'up_pt');
        
    elseif (ceq(condition, {'W', 0, 'new'}))    
        k1h = params.k * cos(params.theta);
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, b, i, j, 1, k1h, 'right_pt', 'up_pt', 'down_pt');
        
    elseif (ceq(condition, {'N', 0, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'N', 0, 'new'}))
        k2h = params.k * sin(params.theta);
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, b, i, j, -1, k2h, 'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'E', 0, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'left_pt', 'up_pt', 'down_pt');
        
    elseif (ceq(condition, {'E', 0, 'new'}))
        k1h = params.k * cos(params.theta);
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, b, i, j, -1, k1h, 'left_pt', 'up_pt', 'down_pt');
    end   
    
end

function [f_central_point] = factory_central_scheme(param)
% function [f_central_point] = factory_central_scheme(param)
% 
% This function is responsible for building the central point scheme; that
% is the point that is not in contact with any boundaries. However this
% point comes with two flavor: the "old" standard 3 points scheme and the
% new scheme (see the publication that is signaled in documentation)
%
% parameters:
%   param: the struct that is describe in the main function
% return:
%   f_central_point: a handle to a function that will compute the
%   coefficient in the scheme matrix.
    % define small function handle that are valids for all the function
    c_std = @(params, A, b, i, j) ctrl_pt_std(params, A, b, i, j);
    c_new = @(params, A, b, i, j) ctrl_pt_new(params, A, b, i, j);     
    
    if strcmp(param.interior, 'new')
        f_central_point = @(params, A, b, i, j) ctrl_pt(params,...
            A, b, i, j, c_new);        
        return
    end
    if strcmp(param.interior,'std')
        f_central_point = @(params, A, b, i, j) ctrl_pt(params,...
            A, b, i, j, c_std);
        return
    end
end

function [A,b] = generic_scheme_creation(params, A, b, i, j, schemes)
%   interior_func,...
%     south_bound, west_bound, north_bound, east_bound,... 
%     sw_corner, nw_corner, ne_corner, se_corner

    m = params.m; 
    n = params.n;
    
%     [i,j]
    
    % south boundary case (without corner)
    if i == 1
        if j==1
            [A,b] = schemes.sw_corner(params, A, b, i, j);
            return
        end
        if j==n
            [A,b] = schemes.se_corner(params, A, b, i, j);
            return
        end
        [A,b] = schemes.south_bound(params, A, b, i, j);
        return
    end
    
    if i == m
       if j==1
           [A,b] = schemes.nw_corner(params, A, b, i, j);
           return
       end
       
       if j==n
            [A,b] = schemes.ne_corner(params, A, b, i, j);
            return
       end
       
       [A,b] = schemes.north_bound(params, A, b, i, j);
       return
    end
    
    if j == 1
        [A,b] = schemes.west_bound(params, A, b, i, j);
        return
    end
    
    if j == m
        [A,b] = schemes.east_bound(params, A, b, i, j);
        return
    end
    
    [A,b] = schemes.interior_func(params, A, b, i, j);
    
end

function [A,b] = left_pt(params, A, b, i)
    [A,b] = generic_pt(params, A, b, i, i-1);
end

function [A,b] = right_pt(params, A, b, i)
    [A,b] = generic_pt(params, A, b, i, i+1);
end

function [A,b] = generic_pt(params, A, b, i, x)
    A(i, x) = A(i, x) - 1;    
end

function [A,b] = ctrl_pt_std(params, A, b, i)
    A(i,i) = 2 - (params.k * params.h).^2;
end

function [A,b] = ctrl_pt_new(params, A, b, i)
    A(i,i) = 2 * cos(params.k * params.h);    
end

function [A,b] = ctrl_pt(params, A, b, i, central_point)
    [A, b] = feval(central_point, params, A, b, i);    
    [A, b] = left_pt(params, A, b, i);
    [A, b] = right_pt(params, A, b, i);        
end

function [A,b] = sommerfeld_side_generic_std(params, A, b, i, j, sgn, int_pt,...
    l_pt, r_pt)
    l = get_scheme_label(params.m, params.n, i, j);    
    kh = params.k * params.h;
    A(l,l) = 4 + sgn * 1i * 2 * kh - (kh).^2;
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(l_pt, params, A, b, i, j);
    [A, b] = feval(r_pt, params, A, b, i, j);
end

function [A,b] = sommerfeld_side_generic_new(params, A, b, i, j, sgn, kxh,...
    int_pt, l_pt, r_pt)
    l = get_scheme_label(params.m, params.n, i, j);    
    kh = params.k * params.h;
    A(l,l) = 4 * besselj(0,kh) + sgn * 2 * 1i * sin(kxh);
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(l_pt, params, A, b, i, j);
    [A, b] = feval(r_pt, params, A, b, i, j);
end

function [A,b] = sommerfeld_generic_corner_new(params, A, b, i, j,...
    sgn1, sgn2, int_pt, ext_pt)
    l = get_scheme_label(params.m, params.n, i, j);
    kh = params.k * params.h;
    k1h = params.k * cos(params.theta);
    k2h = params.k * sin(params.theta);
    A(l,l) = 4 * besselj(0,kh)...
        + sgn1 *  2 * 1i * ( sin(k1h) +  sgn2 * sin(k2h) );
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(ext_pt, params, A, b, i, j);
    [A, b] = feval(ext_pt, params, A, b, i, j);
end

function [A,b] = sommerfeld_generic_corner_std(params, A, b, i, j,...
    int_pt, ext_pt)
    l = get_scheme_label(params.m, params.n, i, j);    
    kh = params.k * params.h;
    A(l,l) = 2 - 1i * sqrt(2) * kh - 0.5 * kh.^2;
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(ext_pt, params, A, b, i, j);
end

function [A,b] = dirichlet_side_generic(params, A, b, i, j,...
    ctrl_pt, dirich_f, int_pt, l_pt, r_pt)
    [A, b] = feval(dirich_f, params, A, b, i, j);    
    [A, b] = feval(ctrl_pt, params, A, b, i, j);
    [A, b] = feval(int_pt, params, A, b, i, j);    
    [A, b] = feval(l_pt, params, A, b, i, j);
    [A, b] = feval(r_pt, params, A, b, i, j);    
end

function [A,b] = dirichlet_corner_generic(params, A, b, i, j,...
    ctrl_pt, dirich_f1, dirich_f2, int_pt1, int_pt2)    
    [A, b] = feval(ctrl_pt, params, A, b, i, j);    
    [A, b] = feval(dirich_f1, params, A, b, i, j);        
    [A, b] = feval(dirich_f2, params, A, b, i, j);    
    [A, b] = feval(int_pt1, params, A, b, i, j);
    [A, b] = feval(int_pt2, params, A, b, i, j);    
end

function [A,b] = dirichlet_wrapper(params, A, b, i, j, func)
    l = get_scheme_label(params.m, params.n, i, j);
    b(l) = b(l) + feval(func, params, A, b, i, j);
end

function check_params(params)
%     classes = {}
%     attributes = {}

    classes = {'struct'};
    attributes = {'nonempty'};
    validateattributes(params, classes, attributes)
        
    classes = {'uint8', 'uint16', 'uint32', 'uint64'};
    attributes = {'noempty', };
    validateattributes(params.m, classes, attributes, 'params.m')
        
    classes = {'numeric'};
    attributes = {'noempty', 'nonnegative'};
    validateattributes(params.k, classes, attributes, 'params.k')
            
    classes = {'numeric'};
    attributes = {'noempty', 'nonnegative'};
    validateattributes(params.k, classes, attributes, 'params.k')
    
    classes = {'char'};
    attributes = {'nonempty'};
    validateattributes(params.interior, classes, attributes, 'params.interior')    
    
    validStrings = {'std', 'new'};
    validatestring(params.interior,validStrings, 'params.interior')   
    
    classes = {'char'};
    attributes = {'nonempty'};
    validateattributes(params.boundary, classes, attributes, 'params.boundary')    
    
    validStrings = {'sommerfeld_new', 'sommerfeld_std', 'dirichlet'};
    validatestring(params.boudary, validStrings, 'params.boundary')   
    
    if ~isfield(params, 'dirichlet')
        error('helmoltz_1d_scheme:argChk',...
            'At least one source must be set (as Dirichlet).')
    end
        
    total_sommerfeld = ...
       ~isfield(params.dirichlet, 'W') && ~isfield(params.dirichlet, 'E');
    if total_sommerfeld
        error('helmoltz_1d_scheme:argChk',...
        'At least one source must be set (as Dirichlet).');                
    end
    
    total_dirichlet = ...
       isfield(params.dirichlet, 'W') && isfield(params.dirichlet, 'E');        
    if ~isfield(params, 'boundary') && ~total_dirichlet
        error('helmoltz_1d_scheme:argChk',...
            'parameter "params.boundary" must be set if not all the border are Dirichlet.')
    end
        
    if ~total_dirichlet && strcmp(params.boundary, 'dirichlet')
        error('helmoltz_two_2d_scheme:argChk',...
        'Either all bounds are Dirichlet either scheme for Sommerfeld computation (''sommerfeld_new'', ''sommerfeld_std'').');                
    end

end


