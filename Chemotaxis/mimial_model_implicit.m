%% Minimal Partial Differenital Equation Model for Chemotaxis%% 

% implicit method
%%
clc
clear 
close all
%%
membrane_size = 5;
grid_size = 100;
num = grid_size^2;
dx = membrane_size/(grid_size - 1); 
x = linspace(0,membrane_size,grid_size);
[dom_x,dom_y]  = meshgrid(x,x);
%%

cx1 = 1.5; cy1 = 2.5;   % Left blob
cx2 = 3.5; cy2 = 2.5;   % Right blob
sigma_u = 0.3;

% Create two Gaussians
u1 = exp(-((dom_x - cx1).^2 + (dom_y - cy1).^2) / (2*sigma_u^2));
u2 = exp(-((dom_x - cx2).^2 + (dom_y - cy2).^2) / (2*sigma_u^2));
u = u1 + u2;

% Initialize v as a smoothed version of u
sigma_v = 1.5;
kernel_size = 7;
[xk, yk] = meshgrid(-kernel_size:kernel_size, -kernel_size:kernel_size);
G = exp(-(xk.^2 + yk.^2) / (2 * sigma_v^2));
G = G / sum(G(:));
v = conv2(u, G, 'same')/10;


u = u(:);
v = v(:);
% [v,~]  = meshgrid(linspace(0,0.1,grid_size),linspace(0,0.1,grid_size));
%%
dt = 0.001;
D = 0.3;
kappa = 0;
%% Operator Consturction
% Need to make sure the boundary conditions is doing its job

[inlap_x, inlap_y] = interial_laplace(grid_size,grid_size);
[div_x, div_y] = interial_divergence(grid_size,grid_size);


neumann_patch = apply_Neumann(sparse(num,num),grid_size,grid_size,'one_sided');
symmetric_patch = apply_Neumann(sparse(num,num),grid_size,grid_size);

lapc = (inlap_x+inlap_y+symmetric_patch)./(2*dx^2);
div_x = (div_x + symmetric_patch)./(2*dx);
div_y = (div_y + symmetric_patch)./(2*dx);

% lapc = apply_Neumann(inlap_x+inlap_y,grid_size,grid_size,'symmetric')./(2*dx^2);
% div_x = apply_Neumann(div_x,grid_size,grid_size,'symmetric')./(2*dx);
% div_y = apply_Neumann(div_y,grid_size,grid_size,'symmetric')./(2*dx);

[upwind_op_x, downwind_op_x] = interial_updownwind(grid_size,grid_size,1);
[upwind_op_y, downwind_op_y] = interial_updownwind(grid_size,grid_size,2);

upwind_op_x = (upwind_op_x+neumann_patch)./dx;
upwind_op_y = (upwind_op_y+neumann_patch)./dx;
downwind_op_x = (downwind_op_x+neumann_patch)./dx;
downwind_op_y = (downwind_op_y+neumann_patch)./dx;

%%
U = zeros(num,1000);
V = zeros(num,1000);


% 
% lapc = laplacian_2D_offdiagonal(grid_size,grid_size,dx);
% lapc = laplacian_neumann_finalize(lapc,grid_size,grid_size,dx);
%%
N = grid_size - 2;


for i = 2:5000



    flux_sign_x = sign(inlap_x*v);
    flux_sign_y = sign(inlap_y*v);
    grad_v_x = (flux_sign_x > 0) .* upwind_op_x*v + (flux_sign_x < 0) .* downwind_op_x*v;
    grad_v_y = (flux_sign_y > 0) .* upwind_op_y*v + (flux_sign_y < 0) .* downwind_op_y*v;

    Jx = kappa * (u.* grad_v_x);
    Jy = kappa * (u.* grad_v_y);

    div_Jv = div_x*Jx+ div_y*Jy;

    % divergence of J

    % flux_sign_x = sign(inlap_x * J);
    % flux_sign_y = sign(inlap_y * J);
    % 
    % selector_x = (flux_sign_x > 0) .* upwind_op_x*J + (flux_sign_x < 0) .* downwind_op_x*J;
    % selector_y = (flux_sign_y > 0) .* upwind_op_y*J + (flux_sign_y < 0) .* downwind_op_y*J;

    % selector_div = selector_x+selector_y;

    u = u + dt*(-div_Jv+ D*lapc*u);

    v = v + dt*(v\lapc+u-v);


    if mod(i,5) == 0
        U(:,i/5) = u;
        V(:,i/5) = v;
    end



end



%%
vid = VideoWriter('UV_concentration_video8 e.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);

for i = 1:1000
    if mod(i,10) == 0
        clf
        subplot(2,1,1)
        z = reshape(U(:,i),grid_size,grid_size);
        imagesc(z); colorbar; axis equal tight;
        title('Concentration of U');

        subplot(2,1,2)
        y = reshape(V(:,i),grid_size,grid_size);
        imagesc(y); colorbar; axis equal tight;
        title('Concentration of V');

        frame = getframe(gcf);
        writeVideo(vid, frame);
        
    end
end 

close(vid);

%%
% vid = VideoWriter('UV_gradient_video.mp4', 'MPEG-4');
% vid.FrameRate = 10;  % Adjust frame rate as needed
% open(vid);
% 
% for i = 1:5000
%     if mod(i,40) == 0
%         clf
%         [dvx, dvy] = gradient(U(:,:,i), 600*dx);
%         quiver(dom_x, dom_y, dvx, dvy,6);
%         axis equal tight;
%         title('Concentration Gradient of U');
% 
% 
%         frame = getframe(gcf);
%         writeVideo(vid, frame);
%     end
% end 
% 
% close(vid);
