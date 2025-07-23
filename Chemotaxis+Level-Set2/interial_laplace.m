function [row_CDM, col_CDM] = interial_laplace(row, col)
% INTERIAL_LAPLACE  Constructs 2D central difference Laplacian operators
% with internal-only (non-boundary) updates in both row and column directions.
%
% Inputs:
%   row - number of rows (y-dimension)
%   col - number of columns (x-dimension)
%
% Outputs:
%   row_CDM - Laplacian operator in row (x) direction using Kronecker product
%   col_CDM - Laplacian operator in column (y) direction using Kronecker product

    %% --- Column-direction Laplacian (∂²/∂y²) ---
    % Construct 1D second-order central difference operator in y-direction
    base_CDM = spdiags([1 -2 1], [-1 0 1], row, row);

    % Zero out top and bottom two rows to exclude boundary updates
    base_CDM(1,:)     = 0;
    base_CDM(2,:)     = 0;
    base_CDM(end,:)   = 0;
    base_CDM(end-1,:) = 0;

    % Identity matrix with zeroed top and bottom rows
    IM_col = speye(col);
    IM_col(1,:)   = 0;
    IM_col(end,:) = 0;

    % Kronecker product: apply base_CDM to each column block
    col_CDM = kron(IM_col, base_CDM);

    %% --- Row-direction Laplacian (∂²/∂x²) ---
    % Construct 1D second-order central difference operator in x-direction
    base_CDM = spdiags([1 -2 1], [-1 0 1], col, col);

    % Zero out left and right two columns to exclude boundary updates
    base_CDM(1,:)     = 0;
    base_CDM(2,:)     = 0;
    base_CDM(end,:)   = 0;
    base_CDM(end-1,:) = 0;

    % Identity matrix with zeroed left and right columns
    IM_row = speye(row);
    IM_row(1,:)   = 0;
    IM_row(end,:) = 0;

    % Kronecker product: apply base_CDM to each row block
    row_CDM = kron(base_CDM, IM_row);
end