function [upwind_op, downwind_op] = exposed_updownwind(row, col, dir)
%EXPOSED_UPDOWNWIND Constructs sparse upwind and downwind operators
%   row, col: grid dimensions
%   dir = 1 → y-direction (rows)
%   dir = 2 → x-direction (columns)
%
%   Both operators are of size (row*col)-by-(row*col), acting on column-stacked fields.

    N = row * col;

    if dir == 2  % y-direction (column-wise traversal)
        % Operators in row direction (vertical difference)
        D_up = spdiags([-1 1], [0 1], row, row);
        % D_up(1, :) = 0;  % open bottom boundary (Neumann)
        % D_up(end, :) = 0;  % open bottom boundary (Neumann)

        D_down = spdiags([-1 1], [-1 0], row, row);
        % D_down(1, :) = 0;  % open top boundary
        % D_down(end, :) = 0;  % open top boundary


        I = speye(col);
        % I([1 end], :) = 0;  % suppress side coupling

        upwind_op   = kron(I, D_up);
        downwind_op = kron(I, D_down);

    elseif dir == 1  % x-direction (row-wise traversal)
        % Operators in col direction (horizontal difference)
        D_up = spdiags([-1 1], [0 1], col, col);

        % D_up(1, :) = 0;  % open right boundary
        % D_up(end, :) = 0;  % open right boundary

        D_down = spdiags([-1 1], [-1 0], col, col);
        % D_down(1, :) = 0;  % open left boundary
        % D_down(end, :) = 0;  % open left boundary


        I = speye(row);
        % No masking needed here since direction is horizontal

        upwind_op   = kron(D_up, I);
        downwind_op = kron(D_down, I);

    else
        error('Invalid direction: dir must be 1 (y) or 2 (x)');
    end

end



