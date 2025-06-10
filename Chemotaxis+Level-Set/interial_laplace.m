function [row_CDM, col_CDM] = interial_laplace(row,col)
    
    col_CDM = spdiags([1 -2 1],[-1 0 1],row,row);
    
    col_CDM(1,:) = 0;
    col_CDM(2,:) = 0;
    col_CDM(end,:) = 0;
    col_CDM(end-1,:) = 0;
    IM = speye(col);
    
    IM(1,:) = 0;
    IM(end,:) = 0;
    
    col_CDM = kron(IM,col_CDM);


    % create central difference for row
    row_CDM = spdiags([1 -2 1],[-1 0 1],col,col);

    row_CDM(1,:) = 0;
    row_CDM(2,:) = 0;
    row_CDM(end,:) = 0;
    row_CDM(end-1,:) = 0;
    IM = speye(row);

    IM(1,:) = 0;
    IM(end,:) = 0;

    row_CDM = kron(row_CDM,IM);

end