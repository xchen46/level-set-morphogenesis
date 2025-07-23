function [row_CDM, col_CDM] = interial_divergence_CDM(row, col)
% INTERIAL_DIVERGENCE_CDM
% Constructs sparse 2D central difference matrices for computing ∂/∂x and ∂/∂y.
% Only internal grid points are updated — boundary rows/columns are zeroed out.
%
% Inputs:
%   row - number of grid points along y-direction (rows)
%   col - number of grid points along x-direction (columns)
%
% Outputs:
%   row_CDM - sparse matrix for ∂/∂x (along rows)
%   col_CDM - sparse matrix for ∂/∂y (along columns)

    num = row * col;  % Total number of grid points

    %% === Y-direction Central Difference Operator (∂/∂y) ===
    % Create 1D central difference stencil in y-direction (along each column)
    base_y = spdiags([-1 0 1], [-1 0 1], row, row);

    % Zero out rows near top and bottom boundaries to avoid invalid stencils
    base_y(1,:)     = 0;
    base_y(2,:)     = 0;
    base_y(end-1,:) = 0;
    base_y(end,:)   = 0;

    % Identity matrix across x-columns, with edges zeroed
    mask_col = speye(col);
    mask_col(1,:)   = 0;
    mask_col(end,:) = 0;

    % Apply operator to each column using Kronecker product
    col_CDM = kron(mask_col, base_y);

    %% === X-direction Central Difference Operator (∂/∂x) ===
    % Create 1D central difference stencil in x-direction (along each row)
    base_x = spdiags([-1 0 1], [-1 0 1], col, col);

    % Zero out rows near left and right boundaries to avoid invalid stencils
    base_x(1,:)     = 0;
    base_x(2,:)     = 0;
    base_x(end-1,:) = 0;
    base_x(end,:)   = 0;

    % Identity matrix across y-rows, with edges zeroed
    mask_row = speye(row);
    mask_row(1,:)   = 0;
    mask_row(end,:) = 0;

    % Apply operator to each row using Kronecker product
    row_CDM = kron(base_x, mask_row);
end
