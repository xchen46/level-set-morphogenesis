function [FW, BW] = interial_divergence_FDM(row, col, dim)
% INTERIAL_DIVERGENCE_FDM
% Constructs sparse forward and backward first-order finite difference operators
% in 2D (∂/∂x or ∂/∂y) with masked boundaries to prevent stencil overflow.
%
% Inputs:
%   row - number of grid points along the y-direction
%   col - number of grid points along the x-direction
%   dim - 1 for ∂/∂x (horizontal), 2 for ∂/∂y (vertical)
%
% Outputs:
%   FW - forward difference operator (size: (row*col) × (row*col))
%   BW - backward difference operator (size: (row*col) × (row*col))
%
% These operators are suitable for use in conservative Laplacian schemes:
%       L = BW * Dr * FW

    if dim == 2  % ∂/∂y (vertical direction)
        % --- Forward difference: f(i+1) - f(i) ---
        base_FW = spdiags([-1 1], [0 1], row, row);
        base_FW([1, end-1:end], :) = 0;  % Remove invalid stencils near top/bottom

        % Mask to disable first and last columns (x-boundary)
        IM_col = speye(col);
        IM_col([1, end], :) = 0;

        % Apply operator down each column
        FW = kron(IM_col, base_FW);

        % --- Backward difference: f(i) - f(i-1) ---
        base_BW = spdiags([-1 1], [-1 0], row, row);
        base_BW([1:2, end-1:end], :) = 0;  % Remove invalid stencils near top/bottom

        % Apply operator down each column
        BW = kron(IM_col, base_BW);

    elseif dim == 1  % ∂/∂x (horizontal direction)
        % --- Forward difference: f(j+1) - f(j) ---
        base_FW = spdiags([-1 1], [0 1], col, col);
        base_FW([1, end-1:end], :) = 0;  % Remove invalid stencils near left/right

        % Mask to disable top and bottom rows (y-boundary)
        IM_row = speye(row);
        IM_row([1, end], :) = 0;

        % Apply operator across each row
        FW = kron(base_FW, IM_row);

        % --- Backward difference: f(j) - f(j-1) ---
        base_BW = spdiags([-1 1], [-1 0], col, col);
        base_BW([1:2, end-1:end], :) = 0;  % Remove invalid stencils near left/right

        % Apply operator across each row
        BW = kron(base_BW, IM_row);
    end
end


