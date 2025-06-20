function int_eye = interial_eye(row,col)
    
    col_CDM = spdiags(1,0,row,row);
    
    col_CDM(1,:) = 0;
    col_CDM(end,:) = 0;
    IM = speye(col);
    
    IM(1,:) = 0;
    IM(end,:) = 0;
    
    int_eye = kron(IM,col_CDM);


    % create central difference for row


end