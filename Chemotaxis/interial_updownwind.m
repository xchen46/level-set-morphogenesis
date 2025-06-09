function [upwind_op, downwind_op] = interial_updownwind(row,col,dir)
    
    if dir == 2 % col
        upwind_op = spdiags([-1 1],[0 1],row,row);

        upwind_op(1,:) = 0;
        upwind_op(2,:) = 0;
        upwind_op(end,:) = 0;
        upwind_op(end-1,:) = 0;
        IM = speye(col);

        IM(1,:) = 0;
        IM(end,:) = 0;

        upwind_op = kron(IM,upwind_op);


        downwind_op = spdiags([-1 1],[-1 0],row,row);

        downwind_op(1,:) = 0;
        downwind_op(2,:) = 0;
        downwind_op(end,:) = 0;
        downwind_op(end-1,:) = 0;

        downwind_op = kron(IM,downwind_op);

    elseif dir == 1 %row
        
        upwind_op = spdiags([-1 1],[0 1],col,col);

        upwind_op(1,:) = 0;
        upwind_op(2,:) = 0;
        upwind_op(end,:) = 0;
        upwind_op(end-1,:) = 0;
        IM = speye(row);

        IM(1,:) = 0;
        IM(end,:) = 0;

        upwind_op = kron(upwind_op,IM);


        downwind_op = spdiags([-1 1],[-1 0],col,col);

        downwind_op(1,:) = 0;
        downwind_op(2,:) = 0;
        downwind_op(end,:) = 0;
        downwind_op(end-1,:) = 0;

        downwind_op = kron(downwind_op,IM);
        % % create central difference for column
        % upwind_op = spdiags([1 -2 1],[-1 0 1],col,col);
        % 
        % upwind_op(1,:) = 0;
        % upwind_op(2,:) = 0;
        % upwind_op(end,:) = 0;
        % upwind_op(end-1,:) = 0;
        % IM = speye(row);
        % 
        % IM(1,:) = 0;
        % IM(end,:) = 0;
        % 
        % row_CDM = kron(upwind_op,IM);
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