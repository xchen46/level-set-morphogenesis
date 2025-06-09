function res = create_mask(M,N)
    mask = ones(M,N);

    % mask(1,:) = 0;
    % mask(end,:) = 0;
    % mask(:,1) = 0;
    % mask(:,end) = 0;
    
    idx_row = [2 2 N-1 N-1];
    idx_col = [2 N-1 2 N-1];
    mask(idx_row,idx_col) = 2;

    idx = 3:N-2;
    mask(2,idx) = 3;
    mask(idx,2) = 3;
    mask(idx,N-1) = 3;
    mask(N-1,idx) = 3;

    [a, b] = ndgrid(idx, idx);
    mask(a,b) = 4;
    
    diag_vals = -mask(:);
    res = spdiags(diag_vals,0, M*N, M*N);
end