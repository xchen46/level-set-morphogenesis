%% Minimal Partial Differenital Equation Model for Level-set function%% 

% in procgress for update
%%
clc
clear 
close all
%% Initialzie Parameter
domain_size = 5;
grid_size = 100;
num = grid_size^2;
dx = domain_size/(grid_size - 1); 
x = linspace(0,domain_size,grid_size);
[phi_x,phi_y]  = meshgrid(x,x);
%%

cx1 = 2.5; cy1 = 2.5;   % Left blob
sigma_u = 0.3;

% Create two Gaussians
phi1 = sqrt((phi_x-cx1).^2+(phi_y-cy1).^2)-1;
velc = 0.1;

% phi1(phi1 < 0) = -1;
% phi1(phi1 > 0) = 1; 
phi = phi1(:);

vx = zeros(grid_size);
vy = zeros(grid_size);
vx(phi_x > 2.5) = phi_x(phi_x > 2.5)-2.5;
vx(phi_x < 2.5) = -phi_x(phi_x < 2.5)+2.5;

vx(phi_x == 2.5) = phi_x(phi_x == 2.5)-2.5;

vy(phi_y > 2.5) = phi_y(phi_y > 2.5)-2.5;
vy(phi_y < 2.5) = -phi_y(phi_y < 2.5)+2.5;

vy(phi_y == 2.5) = phi_y(phi_y == 2.5)-2.5;

vx = vx(:);
vy = vy(:);
%%
dt = 0.001;
v = 0.001;
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

upwind_op_x = (upwind_op_x)./dx;
upwind_op_y = (upwind_op_y)./dx;
downwind_op_x = (downwind_op_x)./dx;
downwind_op_y = (downwind_op_y)./dx;

%%
Heav = zeros(num,1000);
Heav(:,1) = smear_out_heaviside(phi,dx);
% 
% lapc = laplacian_2D_offdiagonal(grid_size,grid_size,dx);
% lapc = laplacian_neumann_finalize(lapc,grid_size,grid_size,dx);
%%
N = grid_size - 2;


for i = 2:5000


    grad_phi_x = (vx < 0) .* upwind_op_x*phi + (vx > 0) .* downwind_op_x*phi;
    grad_phi_y = (vy < 0) .* upwind_op_y*phi + (vy > 0) .* downwind_op_y*phi;

    grad_phi = vx.*grad_phi_x+vy.*grad_phi_y;


    phi = phi - dt*(grad_phi);


    if mod(i,5) == 0
        Heav(:,i/5) = smear_out_heaviside(phi,dx);
    end



end



%%
vid = VideoWriter('UV_concentration_video8.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);

for i = 1:1000
    if mod(i,10) == 0
        clf
        z = reshape(Heav(:,i),grid_size,grid_size);
        contour(z); colorbar; axis equal tight;
        % imagesc(z); colorbar; axis equal tight;
        title('Concentration of U');

        frame = getframe(gcf);
        writeVideo(vid, frame);
        
    end
end 

close(vid);

