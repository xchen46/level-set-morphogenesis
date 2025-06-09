%% Applying Neumann Boundary %%
function S = laplacian_neumann_finalize(S, M, N, dx, dy)

    if nargin < 5, dy = dx; end

    [row, col] = size(S);

    mask = (2/dx^2)*ones(M-2,N-2)+(2/dy^2)*ones(M-2,N-2);


    mask([1 end],:) = mask([1 end],:) - 1/dy^2; % Top and bottom edges
    mask(:,[1 end]) = mask(:, [1 end]) - 1/dx^2;

    mask2 = 1/dx^2*ones(M,N);
    mask2(2:M-1, 2:N-1) = mask;

    diag_vals = -mask2(:);

    idx = sub2ind([row, col], 1:size(S,1), 1:size(S,2));

    S(idx) = diag_vals;

    % Top and bottom boundaries
    for boundary = [1, row - N + 1]
        idx = boundary : boundary + N - 1;
        S(idx, :) = 0; S(:,idx) = 0;
        S(sub2ind([row, col], idx, idx)) = 1;
        S(sub2ind([row, col], idx, idx + (2 * (boundary == 1) - 1) * N)) = -1;
    end

    
end

%%
% function res = apply_neumann(S,M,N)
% 
%     [row, col] = size(S);
% 
%     top = 1:N;
%     bd_top = top+N;
%     btm = row-N+1:row;
%     bd_btm = btm-N;
% 
%     S(top, :) = 0;
%     S(btm, :) = 0;
%     idx = sub2ind([row col], top, bd_top);
%     S(idx) = 1;
%     idx = sub2ind([row col], top, top);
%     S(idx) = -1;
% 
%     idx = sub2ind([row col], btm, bd_btm);
%     S(idx) = 1;
%     idx = sub2ind([row col], btm, btm);
%     S(idx) = -1;
% 
% 
%     inner_1 = N+1:N:row-N;
%     S(inner_1,:) = 0;
%     idx = sub2ind([row col], inner_1, inner_1);
%     S(idx) = -1;
%     idx = sub2ind([row col], inner_1+1, inner_1);
%     S(idx) = S(idx)+1;
% 
%     inner_2 = 2*N:N:row-N;
%     S(inner_2,:) = 0;
%     idx = sub2ind([row col], inner_2, inner_2);
%     S(idx) = -1;
%     idx = sub2ind([row col], inner_2-1, inner_2);
%     S(idx) = S(idx)+1;
% 
%     res = S;
% 
% end


% % Left and right boundaries (excluding top/bottom)
    % for offset = [1,2]
    %     if offset == 1
    %         inner = N+1:N:row-N;
    %         S(sub2ind([row, col], inner, inner)) = -1;
    %         S(sub2ind([row, col], inner+1, inner)) = 2;
    %     else
    %         inner = 2*N:N:row-N;
    %         S(sub2ind([row, col], inner, inner)) = -1;
    %         S(sub2ind([row, col], inner-1, inner)) = 2;
    %     end

        % S(sub2ind([row, col], inner + (2 * (offset == 1) - 1))* N,inner) = 2;
    % end