%% Kronecker sum of discrete Laplacians %%

% this creates 2D laplacian with no diagonal
% and row boundary ready to be enforced
% without applying the boundary condition

function [res] = laplacian_2D_offdiagonal(M,N,dx,dy)
    % If dy is not provided, assume dy = dx
    if nargin < 4, dy = dx; end

    Dx = spdiags([1 2 1],[-1 0 1],N,N);
    % Dx(2,1) = 0; Dx(N-1,N) = 0;
    Dy = spdiags([1 2 1],[-1 0 1],M,M); 
    
    IN = speye(N); 
    IN(:, [1 N]) = 0;
    IM = speye(M);
    
    res = (1/dx^2)*kron(IM,Dx)+ (1/dy^2)*kron(Dy,IN);
    
end