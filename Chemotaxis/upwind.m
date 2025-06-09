function res = upwind(U,dx,i,j,up,idx)
    
    
    if idx == 0
        if up == 1
            res = (U(i,j)-U(i-1,j))/dx;
        elseif up == 0
            res = (U(i+1,j)-U(i,j))/dx;
        end
    elseif idx == 1
        if up == 1
            res = (U(i,j)-U(i,j-1))/dx;
        elseif up == 0
            res = (U(i,j+1)-U(i,j))/dx;
        end
    else
        error('upwind: invalid dimension index (should be 0 or 1)');
    end


end