function res = CDM(U,dx,i,j,idx)
    if idx == 0
        res = (U(i+1,j)+U(i-1,j)-2*U(i,j))/dx^2;
    elseif idx == 1
        res = (U(i,j+1)+U(i,j-1)-2*U(i,j))/dx^2;
    else
        print('CDM, not correct dimension index')
    end
end