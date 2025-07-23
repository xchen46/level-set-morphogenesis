function [row_CDM, col_CDM] = full_divergence(row,col,dx,mode)

% CENTRAL_DIFFERENCE_2D Generate 2D central difference matrices for ∂/∂x and ∂/∂y
%
%   [Dx, Dy] = central_difference_2D(row, col, dx, mode)
%
%   Inputs:
%       row, col - number of grid points in y and x
%       dx       - grid spacing
%       mode     - 'bound' or 'unbound' (Neumann patch)
%
%   Outputs:
%       Dx - sparse matrix for ∂/∂x
%       Dy - sparse matrix for ∂/∂y

    if nargin < 4
        mode = 'unbound';
    end

    num = row * col;
    
    % Y-derivative (∂/∂y)
    col_CDM = spdiags([-1 0 1],[-1 0 1],row,row);
    
    % col_CDM(1,:) = 0;
    % col_CDM(2,:) = 0;
    % col_CDM(end,:) = 0;
    % col_CDM(end-1,:) = 0;
    IM = speye(col);
    
    % Zero out boundary rows to disable invalid stencil application near boundaries
    % IM(1,:) = 0; IM(end,:) = 0;
    
    col_CDM = kron(IM,col_CDM);


    % X-derivative (∂/∂x)
    row_CDM = spdiags([-1 0 1],[-1 0 1],col,col);

    % row_CDM(1,:) = 0;
    % row_CDM(2,:) = 0;
    % row_CDM(end,:) = 0;
    % row_CDM(end-1,:) = 0;

    % Zero out boundary rows to disable invalid stencil application near boundaries
    IM = speye(row);
    % IM(1,:) = 0;
    % IM(end,:) = 0;

    row_CDM = kron(row_CDM,IM);

    if strcmp(mode, 'bound')

        symmetric_patch = apply_Neumann(sparse(num,num),row,col);
        row_CDM = (row_CDM+symmetric_patch)./(2*dx);
        col_CDM = (col_CDM+symmetric_patch)./(2*dx);
    elseif strcmp(mode, 'unbound')
        row_CDM = (row_CDM)./(2*dx);
        col_CDM = (col_CDM)./(2*dx);
    end

end