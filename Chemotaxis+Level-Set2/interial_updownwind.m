function [upwind_op, downwind_op] = interial_updownwind(row, col, dir)
% INTERIAL_UPDOWNWIND
% Constructs sparse upwind and downwind operators (first-order FDM) for a 2D grid.
% Operators are masked to act only on interior points.
%
% Inputs:
%   row - number of grid points in y-direction
%   col - number of grid points in x-direction
%   dir - direction: 1 for x (horizontal), 2 for y (vertical)
%
% Outputs:
%   upwind_op   - forward-biased difference operator (e.g., f(i+1) - f(i))
%   downwind_op - backward-biased difference operator (e.g., f(i) - f(i-1))

    if dir == 2  % === Vertical direction (∂/∂y) ===
        % --- Upwind operator: f(i+1) - f(i) ---
        base_up = spdiags([-1 1], [0 1], row, row);
        base_up([1, end-1:end], :) = 0;  % Mask top and bottom rows

        % Identity over columns (horizontal)
        mask_col = speye(col);
        mask_col([1:2, end-1:end], :) = 0;  % Mask near side boundaries

        upwind_op = kron(mask_col, base_up);

        % --- Downwind operator: f(i) - f(i-1) ---
        base_down = spdiags([-1 1], [-1 0], row, row);
        base_down([1:2, end], :) = 0;  % Mask invalid backward rows

        downwind_op = kron(mask_col, base_down);

    elseif dir == 1  % === Horizontal direction (∂/∂x) ===
        % --- Upwind operator: f(j+1) - f(j) ---
        base_up = spdiags([-1 1], [0 1], col, col);
        base_up([1:2, end-1:end], :) = 0;  % Mask near boundaries

        % Identity over rows (vertical)
        mask_row = speye(row);
        mask_row([1, end], :) = 0;

        upwind_op = kron(base_up, mask_row);

        % --- Downwind operator: f(j) - f(j-1) ---
        base_down = spdiags([-1 1], [-1 0], col, col);
        base_down([1:2, end-1:end], :) = 0;  % Mask near boundaries

        downwind_op = kron(base_down, mask_row);
    end
end




% function [upwind, downwind] = updownwind(primary_dim,orthogonal_dim,d)
%     % spacing is column for x
%     % spacing is row for y
% 
%     % need to fix the y directions!!!!!!! !!!!!!
%     up = (1/d)*spdiags([1 -1],[0 1],orthogonal_dim,orthogonal_dim);
%     % up(end,end-1) = -1;
% 
%     up(1,:) = 0;
%     up(end,:) = 0;
% 
%     IM = speye(primary_dim); 
% 
%     IM(1,:) = 0;
%     IM(end,:) = 0;
% 
%     upwind = kron(IM,up);
% 
%     % down = (1/d)*spdiags([1 -1],[-1 0],primary_dim,primary_dim);
%     % % down(1,2) = -1;
%     % 
%     % down(1,:) = 0;
%     % down(end,:) = 0;
%     % downwind = kron(IM,down);
% 
% 
% 
%     % downwind = (1/dy)*spdiags([1 -1],[0 1],N,N);
% 
%     % [row, ~] = size(upwind);
%     % 
%     % % Top and bottom boundaries
%     % for boundary = [1, row - primary_dim + 1]
%     %     idx = boundary : boundary + primary_dim - 1;
%     %     upwind(idx, :) = 0; upwind(:,idx) = 0;
%     %     % upwind(sub2ind([row, col], idx, idx)) = -1/d^2;
%     %     % upwind(sub2ind([row, col], idx, idx + (2 * (boundary == 1) - 1) * primary_dim)) = 1/d^2;
%     %     % 
%     %     % downwind(idx, :) = 0; downwind(:,idx) = 0;
%     %     % downwind(sub2ind([row, col], idx, idx)) = -1/d^2;
%     %     % downwind(sub2ind([row, col], idx, idx + (2 * (boundary == 1) - 1) * primary_dim)) = 1/d^2;
%     % end
% 
% end