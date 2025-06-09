%% easy laplace model for sparse matrix play

%%

membrane_size = 5;
grid_size = 200;
x = linspace(0,membrane_size,grid_size);
dx = membrane_size/(grid_size - 1); 
[dom_x,dom_y]  = meshgrid(x,x);
%%
D = 0.1;
kappa = 1;
u = bolb_rand(x,membrane_size);
u = u(:);

U = zeros(length(u),5000);
U(:,1) = u;
dt = 0.0001;
%% construct the laplace operator
lapc = laplacian_2D_offdiagonal(grid_size,grid_size,dx);
lapc = laplacian_neumann_finalize(lapc,grid_size,grid_size,dx);

for i = 2:5000
    u = U(:,i-1);
    
    w = u + D*dt*lapc*u;
    
    U(:,i) = w;

end

%%
for i = 1:5000
    if mod(i,40) == 0
        clf
        u = reshape(U(:,i),grid_size,grid_size);
        imagesc(u); colorbar; axis equal tight;

        pause(0.1);
    end
end 