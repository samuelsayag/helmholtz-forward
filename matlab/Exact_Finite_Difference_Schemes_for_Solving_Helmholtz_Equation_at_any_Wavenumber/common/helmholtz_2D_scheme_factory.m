function [ func_scheme ] = helmholtz_2D_scheme_factory( params )
%SCHEME_FACTORY 
% The purpose of this function is to offer a flexible way of creating a
% scheme of the Helmholtz equation. The scheme is computed by considering a
% rectangular region. One side is a source (the computation of the matrix
% will be similar as a Dirichlet Boundary Condition) but the function of
% the source is left to the client. The three other sides may be
% constrained by Dirichlet Boundary Condition or Sommerfeld condition.
%
% return: 
%   func_scheme: a function that may directly be used as a scheme to create 
% the matrix.n
% parameters:
%   param: a struct that contains
%       m: number of line 
%       n: number of column
%       k: 
%       h: step choosen for the scheme
%       theta: angle (rd) to compute Sommerfeld boundary
%       interior: 'new', 'std'
%       boundary: 'new', 'std'
%       dirichlet.W: (optional) function handle
%       dirichlet.N: (optional) function handle   
%       dirichlet.E: (optional) function handle
%       dirichlet.S: (optional) function handle
%       Note:
%           - a Dirichlet function is of signature: f(A, b, i, j) with 'A'
% the scheme matrix ((m*n)²), 'b' the vector such that the solution x of 
% A x = b is the solution of the Helmholtz equation.
%           - at least one of the side must be of Dirichlet Boundary kind.
%           - if a Dirichlet function is not given the side is considered 
% under a Sommerfeld Boundary.

% test if the south_source parameter is a well defined function
check_params(params)

% preserve not implemented mechanism (for the moment)
params.n = params.n;
% dummy value to sommerfeld boundary in case all is dirichlet. The value
% is tested later and may cause "null pointer".
if ~isfield(params , 'boundary')
    params.boundary = 'dirichlet';
end

dirich_s = isa(params.dirichlet.S, 'function_handle'); 
dirich_w = isa(params.dirichlet.W, 'function_handle'); 
dirich_n = isa(params.dirichlet.N, 'function_handle'); 
dirich_e = isa(params.dirichlet.E, 'function_handle'); 

% THROW EXCEPTION IF ALL ARE FALSE (impossible to apply Sommerfeld 
% on all the sides

schemes.interior_func = factory_central_scheme(params);

key = {'S', dirich_s, params.interior, params.boundary};
schemes.south_bound = factory_side_scheme(key, params);

key = {'W', dirich_w, params.interior, params.boundary};
schemes.west_bound =  factory_side_scheme(key, params);

key = {'N', dirich_n, params.interior, params.boundary};
schemes.north_bound = factory_side_scheme(key, params);

key = {'E', dirich_e, params.interior, params.boundary};
schemes.east_bound =  factory_side_scheme(key, params);

key = {'N', dirich_n, 'E', dirich_e, params.interior, params.boundary};
schemes.ne_corner = factory_corner_scheme(key, params);

key = {'E', dirich_e, 'S', dirich_s, params.interior, params.boundary};
schemes.se_corner = factory_corner_scheme(key, params);

key = {'S', dirich_s, 'W', dirich_w, params.interior, params.boundary};
schemes.sw_corner = factory_corner_scheme(key, params);

key = {'W', dirich_w, 'N', dirich_n, params.interior, params.boundary};
schemes.nw_corner = factory_corner_scheme(key, params);

func_scheme = @(params, A, b, i, j) generic_scheme_creation(...
    params, A, b, i, j, schemes);

end

function [f_corn_pt] = factory_corner_scheme(key, params)
% sides are given in the right direction 
    % define small function handle that are valids for all the function
    c_std = @(params, A, b, i, j)ctrl_pt_std(params, A, b, i, j);
    c_new = @(params, A, b, i, j)ctrl_pt_new(params, A, b, i, j);    
    % define an equal operator for dirichlet choice    
    ceq = @(x,y) strcmp(x{1},y{1}) && x{2} == y{2} && strcmp(x{3},y{3})...
        && x{4} == y{4} && strcmp(x{5},y{5});
    
    
    condition = key(1:5);
   % dirichlet on both sides
   if(ceq(condition, {'N', 1, 'E', 1, 'std'}))
       f_corn_pt = @(params, A, b, i, j) dirichlet_corner_generic(...
           params, A, b, i, j, c_std, params.dirichlet.N,...
           params.dirichlet.E, 'left_pt', 'down_pt');
       
   elseif (ceq(condition, {'N', 1, 'E', 1, 'new'}))
       f_corn_pt = @(params, A, b, i, j) dirichlet_corner_generic(...
           params, A, b, i, j, c_new, params.dirichlet.N,...
           params.dirichlet.E, 'left_pt', 'down_pt');
       
   elseif (ceq(condition, {'E', 1, 'S', 1, 'std'}))
       f_corn_pt = @(params, A, b, i, j)dirichlet_corner_generic(...
           params, A, b, i, j, c_std, params.dirichlet.S,...
           params.dirichlet.E, 'left_pt', 'up_pt');
       
   elseif (ceq(condition, {'E', 1, 'S', 1, 'new'}))
       f_corn_pt = @(params, A, b, i, j)dirichlet_corner_generic(...
           params, A, b, i, j, c_new, params.dirichlet.S,...
           params.dirichlet.E, 'left_pt', 'up_pt');
       
   elseif (ceq(condition, {'S', 1, 'W', 1, 'std'}))
       f_corn_pt = @(params, A, b, i, j) dirichlet_corner_generic(...
           params, A, b, i, j, c_std, params.dirichlet.S,...
           params.dirichlet.W, 'right_pt', 'up_pt');
       
   elseif (ceq(condition, {'S', 1, 'W', 1, 'new'}))
       f_corn_pt = @(params, A, b, i, j) dirichlet_corner_generic(...
           params, A, b, i, j, c_new, params.dirichlet.S,...
           params.dirichlet.W, 'right_pt', 'up_pt');
       
   elseif (ceq(condition, {'W', 1, 'N', 1, 'std'}))
       f_corn_pt = @(params, A, b, i, j) dirichlet_corner_generic(...
           params, A, b, i, j, c_std, params.dirichlet.W,...
           params.dirichlet.N, 'right_pt', 'down_pt');
       
   elseif (ceq(condition, {'W', 1, 'N', 1, 'new'}))
       f_corn_pt = @(params, A, b, i, j) dirichlet_corner_generic(...
           params, A, b, i, j, c_new, params.dirichlet.W,...
           params.dirichlet.N, 'right_pt', 'down_pt');
   
   end
    
    % define an equal operator for sommerfeld choice
    condition = key([1,2,3,4,6]) ;
   
    % sommerfeld on both sides
    if(ceq(condition, {'N', 1, 'E', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_std(...
            params, A, b, i, j, 'down_pt','left_pt');
        
    elseif (ceq(condition, {'N', 1, 'E', 1, 'new'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_new(...
            params, A, b, i, j, -1, 1, 'down_pt', 'left_pt');
        
    elseif (ceq(condition, {'E', 1, 'S', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_std(...
            params, A, b, i, j, 'up_pt','left_pt');
        
    elseif (ceq(condition, {'E', 1, 'S', 1, 'new'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_new(...
            params, A, b, i, j, -1, -1, 'up_pt', 'left_pt');
        
    elseif (ceq(condition, {'S', 1, 'W', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_std(...
            params, A, b, i, j, 'up_pt','right_pt');
        
    elseif (ceq(condition, {'S', 1, 'W', 1, 'new'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_new(...
            params, A, b, i, j, 1, 1, 'up_pt', 'right_pt');
        
    elseif (ceq(condition, {'W', 1, 'N', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_std(...
            params, A, b, i, j, 'down_pt','right_pt');
        
    elseif (ceq(condition, {'W', 1, 'N', 1, 'new'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_generic_corner_new(...
            params, A, b, i, j, 1, -1, 'down_pt', 'right_pt');
      
    end
  
   
    % define an equal operator for DIRICHLET-SOMMERFELD (repectively) choice
    condition = key([1,2,3,4,6]);
    
    if(ceq(condition, {'N', 1, 'E', 0, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'left_pt', params.dirichlet.N, 'down_pt');
        
    elseif (ceq(condition, {'N', 1, 'E', 0, 'new'}))
        k1h = params.k * cos(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, -1, k1h, 'left_pt', params.dirichlet.N, 'down_pt');
        
    elseif (ceq(condition, {'E', 1, 'S', 0, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, 1, 'up_pt', params.dirichlet.E, 'left_pt');
        
    elseif (ceq(condition, {'E', 1, 'S', 0, 'new'}))
        k2h = params.k * sin(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, 1, k2h, 'up_pt', params.dirichlet.E, 'left_pt');
        
    elseif (ceq(condition, {'S', 1, 'W', 0, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, 1, 'right_pt', params.dirichlet.S, 'up_pt');
        
    elseif (ceq(condition, {'S', 1, 'W', 0, 'new'}))
        k1h = params.k * cos(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, 1, k1h, 'right_pt', params.dirichlet.S, 'up_pt');
        
    elseif (ceq(condition, {'W', 1, 'N', 0, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'down_pt', params.dirichlet.W, 'right_pt');
        
    elseif (ceq(condition, {'W', 1, 'N', 0, 'new'}))
        k2h = params.k * sin(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, -1, k2h, 'down_pt', params.dirichlet.W, 'right_pt');
    end
   
    % define an equal operator for SOMMERFELD-DIRICHLET choice
    condition = key([1,2,3,4,6]);
    
    if(ceq(condition, {'N', 0, 'E', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'down_pt', params.dirichlet.E, 'left_pt');
        
    elseif (ceq(condition, {'N', 0, 'E', 1, 'new'}))
        k2h = params.k * sin(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, -1, k2h, 'down_pt', params.dirichlet.E, 'left_pt');
        
    elseif (ceq(condition, {'E', 0, 'S', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'left_pt', params.dirichlet.S, 'up_pt');
        
    elseif (ceq(condition, {'E', 0, 'S', 1, 'new'}))
        k1h = params.k * cos(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, -1, k1h, 'left_pt', params.dirichlet.S, 'up_pt');
        
    elseif (ceq(condition, {'S', 0, 'W', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, 1, 'up_pt', params.dirichlet.W, 'right_pt');
        
    elseif (ceq(condition, {'S', 0, 'W', 1, 'new'}))
        k2h = params.k * sin(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, 1, k2h, 'up_pt', params.dirichlet.W, 'right_pt');
        
    elseif (ceq(condition, {'W', 0, 'N', 1, 'std'}))
        f_corn_pt = @(params, A, b, i, j) sommerfeld_side_generic_std(...
            params, A, b, i, j, 1, 'right_pt', params.dirichlet.N, 'down_pt');
        
    elseif (ceq(condition, {'W', 0, 'N', 1, 'new'}))
        k1h = params.k * cos(params.theta);
        f_corn_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, 1, k1h, 'right_pt', params.dirichlet.N, 'down_pt');
    end
    
end

function [f_ctrl_pt] = factory_side_scheme(key, params)
    % define small function handle that are valids for all the function
    c_std = @(params, A, b, i, j)ctrl_pt_std(params, A, b, i, j);
    c_new = @(params, A, b, i, j)ctrl_pt_new(params, A, b, i, j);
    %dirichlet    
    % define an equal operator for dirichlet choice
    condition = key(1:3);
    ceq = @(x,y) strcmp(x{1}, y{1}) && (x{2} == y{2}) && strcmp(x{3},y{3});
    if (ceq(condition, {'S', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j) dirichlet_side_generic(...
            params, A, b, i, j, c_std, params.dirichlet.S,...
            'up_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'S', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, params.dirichlet.S,...
            'up_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'W', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_std, params.dirichlet.W,...
            'down_pt', 'up_pt', 'right_pt');
        
    elseif (ceq(condition, {'W', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, params.dirichlet.W,...
            'down_pt', 'up_pt', 'right_pt');
        
    elseif (ceq(condition, {'N', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_std, params.dirichlet.N,...
            'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'N', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, params.dirichlet.N,...
            'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'E', 1, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_std, params.dirichlet.E,...
            'down_pt', 'left_pt', 'up_pt');
        
    elseif (ceq(condition, {'E', 1, 'new'}))
        f_ctrl_pt = @(params, A, b, i, j)dirichlet_side_generic(...
            params, A, b, i, j, c_new, params.dirichlet.E,...
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
            params, A, i, j, 1, k2h, 'up_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'W', 0, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_std(...
            params, A, b, i, j, 1, 'right_pt', 'down_pt', 'up_pt');
        
    elseif (ceq(condition, {'W', 0, 'new'}))    
        k1h = params.k * cos(params.theta);
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, 1, k1h, 'right_pt', 'up_pt', 'down_pt');
        
    elseif (ceq(condition, {'N', 0, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'N', 0, 'new'}))
        k2h = params.k * sin(params.theta);
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, -1, k2h, 'down_pt', 'left_pt', 'right_pt');
        
    elseif (ceq(condition, {'E', 0, 'std'}))
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_std(...
            params, A, b, i, j, -1, 'left_pt', 'up_pt', 'down_pt');
        
    elseif (ceq(condition, {'E', 0, 'new'}))
        k1h = params.k * cos(params.theta);
        f_ctrl_pt = @(params, A, b, i, j)sommerfeld_side_generic_new(...
            params, A, i, j, -1, k1h, 'left_pt', 'up_pt', 'down_pt');
    end   
    
end

function [f_central_point] = factory_central_scheme(param)
% function [f_central_point] = factory_central_scheme(param)
% 
% This function is responsible for building the central point scheme; that
% is the point that is not in contact with any boundaries. However this
% point comes with two flavor: the "old" standard 5 points scheme and the
% new scheme (see the publication that is signaled in documentation)
%
% parameters:
%   param: the struct that is describe in the main function
% return:
%   f_central_point: a handle to a function that will handle the
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
    
    [i,j]
    
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

function [label] = get_scheme_label(m, n, i, j)
    label = j + (m - i) * n;
end

function [A,b] = left_pt(params, A, b, i, j)
    [A,b] = generic_pt(params, A, b, i, j, i, j-1);
end

function [A,b] = right_pt(params, A, b, i, j)
    [A,b] = generic_pt(params, A, b, i, j, i, j+1);
end

function [A,b] = up_pt(params, A, b, i, j)
    [A,b] = generic_pt(params, A, b, i, j, i+1, j);
end

function [A,b] = down_pt(params, A, b, i, j)
    [A,b] = generic_pt(params, A, b, i, j, i-1, j);
end

function [A,b] = generic_pt(params, A, b, i, j, x, y)
    l = get_scheme_label(params.m, params.n, i, j);
    c = get_scheme_label(params.m, params.n, x, y);
    A(l, c) = A(l, c) - 1;    
end

function [A,b] = ctrl_pt_std(params, A, b, i, j)
    l = get_scheme_label(params.m, params.n, i, j);    
    A(l,l) = 4 - (params.k * params.h).^2;
end

function [A,b] = ctrl_pt_new(params, A, b, i, j)
    l = get_scheme_label(params.m, params.n, i, j);    
    A(l,l) = 4 * besselj(0, params.k * params.h);    
end

function [A,b] = ctrl_pt(params, A, b, i, j, central_point)
    [A, b] = feval(central_point, params, A, b, i, j);    
    [A, b] = left_pt(params, A, b, i, j);
    [A, b] = up_pt(params, A, b, i, j);    
    [A, b] = right_pt(params, A, b, i, j);    
    [A, b] = down_pt(params, A, b, i, j);    
end

function sommerfeld_side_generic_std(params, A, b, i, j, sgn, int_pt,...
    l_pt, r_pt)
    l = get_scheme_label(params.m, params.n, i, j);    
    kh = params.k * params.h;
    A(l,l) = 4 + sgn * 1i * 2 * kh - (kh).^2;
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(l_pt, params, A, b, i, j);
    [A, b] = feval(r_pt, params, A, b, i, j);
end

function [ A ] = sommerfeld_side_generic_new(params, A, b, i, j, sgn, kxh,...
    int_pt, l_pt, r_pt)
    l = get_scheme_label(params.m, params.n, i, j);    
    kh = params.k * params.h;
    A(l,l) = 4 * besselj(0,kh) + sgn * 2 * 1i * sin(kxh);
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(l_pt, params, A, b, i, j);
    [A, b] = feval(r_pt, params, A, b, i, j);
end

function [A] = sommerfeld_generic_corner_new(params, A, b, i, j,...
    sgn1, sgn2, int_pt, ext_pt)
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

function [A] = sommerfeld_generic_corner_std(params, A, b, i, j,...
    int_pt, ext_pt)
    kh = params.k * params.h;
    A(l,l) = 2 - 1i * sqrt(2) * kh - 0.5 * kh.^2;
    [A, b] = feval(int_pt, params, A, b, i, j);
    [A, b] = feval(ext_pt, params, A, b, i, j);
end

function [A, b] = dirichlet_side_generic(params, A, b, i, j,...
    ctrl_pt, dirich_f, int_pt, l_pt, r_pt)
    
    dirich = @(params,A,b,i,j) dirichlet_wrapper(params,...
        A, b, i, j, dirich_f);

    [A, b] = feval(dirich, params, A, b, i, j);    
    [A, b] = feval(ctrl_pt, params, A, b, i, j);
    [A, b] = feval(int_pt, params, A, b, i, j);    
    [A, b] = feval(l_pt, params, A, b, i, j);
    [A, b] = feval(r_pt, params, A, b, i, j);    
end

function [A, b] = dirichlet_corner_generic(params, A, b, i, j,...
    ctrl_pt, dirich_f1, dirich_f2, int_pt1, int_pt2)    
    
    dirich1 = @(params,A,b,i,j) dirichlet_wrapper(params,...
        A, b, i, j, dirich_f1);
    dirich2 = @(params,A,b,i,j) dirichlet_wrapper(params,...
        A, b, i, j, dirich_f2);
    
    [A, b] = feval(ctrl_pt, params, A, b, i, j);    
    [A, b] = feval(dirich1, params, A, b, i, j);        
    [A, b] = feval(dirich2, params, A, b, i, j);    
    [A, b] = feval(int_pt1, params, A, b, i, j);
    [A, b] = feval(int_pt2, params, A, b, i, j);    
end

function [A,b] = dirichlet_wrapper(params, A, b, i, j, func)
    l = get_scheme_label(params.m, params.n, i, j);
    b(l) = b(l) + feval(func, params, A, b, i, j);
end

function check_params(params)
    if nargin ~= 1
        error('helmoltz_two_2d_scheme:argChk',...
            'Need a struct as param (see documentation).')
    end
    
    if ~isa(params, 'struct')
        error('helmoltz_two_2d_scheme:argChk',...
            'the parameter "params" of the main function is not a struct.')
    end
    
    if ~isa(params, 'struct')
        error('helmoltz_two_2d_scheme:argChk',...
            'the parameter "params" of the main function is not a struct.')
    end
    
    if ~isfield(params, 'm')
        error('helmoltz_two_2d_scheme:argChk',...
            'parameter "params.m" must be set and is obligatory.')
    end
    
    if ~isfield(params, 'n')
        error('helmoltz_two_2d_scheme:argChk',...
            'parameter "params.n" must be set and is obligatory.')
    end
    
    if ~isfield(params, 'k')
        error('helmoltz_two_2d_scheme:argChk',...
            'parameter "params.k" must be set and is obligatory.')
    end
    
    if ~isfield(params, 'h')
        error('helmoltz_two_2d_scheme:argChk',...
            'parameter "params.h" must be set and is obligatory.')
    end
    
    if ~isfield(params, 'interior')
        error('helmoltz_two_2d_scheme:argChk',...
            'parameter "params.interior" must be set and is obligatory.')
    end
    
    if ~isfield(params, 'dirichlet')
        error('helmoltz_two_2d_scheme:argChk',...
            'At least one source must be set (as Dirichlet).')
    end
        
    total_sommerfeld = ...
       ~isfield(params.dirichlet, 'N') && ~isfield(params.dirichlet, 'E')...
    && ~isfield(params.dirichlet, 'S') && ~isfield(params.dirichlet, 'W');
    if total_sommerfeld
        error('helmoltz_two_2d_scheme:argChk',...
        'At least one source must be set (as Dirichlet).');                
    end
    
    total_dirichlet = ...
       isfield(params.dirichlet, 'N') && isfield(params.dirichlet, 'E')...
    && isfield(params.dirichlet, 'S') && isfield(params.dirichlet, 'W');    
    
    if ~isfield(params, 'boundary') && ~total_dirichlet
        error('helmoltz_two_2d_scheme:argChk',...
            'parameter "params.boundary" must be set if not all the border are Dirichlet.')
    end

    scheme_type_cond = strcmp(params.interior, 'new') || ...
      strcmp(params.interior, 'std');
    if ~scheme_type_cond
        error('helmoltz_two_2d_scheme:argChk',...
        'params.interior must be set to "std" or "new"');        
    end

    if ~isfield(params, 'theta') && ~isfield(params, 'dirichlet')...
            && (strcmp(params.boundary, 'new') || ...
            strcmp(params.boundary, 'std')) 
        error('helmoltz_two_2d_scheme:argChk',...
        'Either all bounds are Dirichlet either specify theta for Sommerfeld computation.');                
    end    
    
    total_dirichlet =  ...
        isfield(params.dirichlet, 'N') && isfield(params.dirichlet, 'E')...
     && isfield(params.dirichlet, 'S') && isfield(params.dirichlet, 'W');
    
    if ~isfield(params, 'theta') && ~total_dirichlet...
            && strcmp(params.boundary, 'new')
        error('helmoltz_two_2d_scheme:argChk',...
        'Either all bounds are Dirichlet either specify theta for Sommerfeld computation.');                
    end
    

 
end

