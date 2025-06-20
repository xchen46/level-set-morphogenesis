%% Minimal Partial Differenital Equation Model for Level-Set + Chemotaxis%% 

% in progress for update
%% clear all the workspace, close all the windows, and clear all command
clc
clear 
close all
%% Initialize Parameters
membrane_size = 10; % Membrane Size
grid_size = 200; % Number of Grid_Point
num = grid_size^2; % Total Number of entries
dx = membrane_size/(grid_size - 1);  % Grid step
x = linspace(0,membrane_size,grid_size); 
[dom_x,dom_y]  = meshgrid(x,x); % Create grid

%% Operator Consturction
% Need to make sure the boundary conditions is doing its job

neumann_patch = apply_Neumann(sparse(num,num),grid_size,grid_size,'one_sided');
symmetric_patch = apply_Neumann(sparse(num,num),grid_size,grid_size,'symmetric');

[inlap_x, inlap_y] = interial_laplace(grid_size,grid_size);
lapc = (inlap_x+inlap_y+symmetric_patch)./(2*dx^2);

[Dx, Dy] = interial_divergence(grid_size,grid_size,dx,'unbound');
Dx = Dx + neumann_patch./dx;
Dy = Dy + neumann_patch./dx;

[upwind_op_x, downwind_op_x] = interial_updownwind(grid_size,grid_size,1);
[upwind_op_y, downwind_op_y] = interial_updownwind(grid_size,grid_size,2);

upwind_op_x = (upwind_op_x+neumann_patch)./dx;
upwind_op_y = (upwind_op_y+neumann_patch)./dx;
downwind_op_x = (downwind_op_x+neumann_patch)./dx;
downwind_op_y = (downwind_op_y+neumann_patch)./dx;

%%

[u,v] = initalize_two_blob(dom_x,dom_y,4,6,1.5,1.5); % initialize u and v

phi = dom_y-2.5;
phi = phi(:);


%% Simulation Parameters

dt =  0.001;
D = 0.2;
kappa = 0.3;
velc = 0.2;

U = zeros(num,1000);
V = zeros(num,1000);

Rho = zeros(num,1000);
Rho(:,1) = 1-smear_out_heaviside(phi,dx);

u = Rho(:,1).*u;
v = Rho(:,1).*v;

%%

for i = 2:5000
    
    % Simulating what is inside the boundary
    flux_sign_x = sign(inlap_x*v);
    flux_sign_y = sign(inlap_y*v);
    grad_v_x = (flux_sign_x > 0) .* upwind_op_x*v + (flux_sign_x < 0) .* downwind_op_x*v;
    grad_v_y = (flux_sign_y > 0) .* upwind_op_y*v + (flux_sign_y < 0) .* downwind_op_y*v;

    Jx = kappa * (u.* grad_v_x);
    Jy = kappa * (u.* grad_v_y);

    div_Jv = Dx*Jx+ Dy*Jy;


    u = u + dt*(-div_Jv+ D*lapc*u);

    v = v + dt*(lapc*v+u-v);


    if mod(i,5) == 0
        U(:,i/5) = u;
        V(:,i/5) = v;
        Rho(:,i/5) = rho;
    end



    % Simulating rho
    [nx, ny] = phi_normal(phi,Dx,Dy);
    dphi_dx = (nx < 0) .* upwind_op_x*phi + (nx > 0) .* downwind_op_x*phi;
    dphi_dy = (ny < 0) .* upwind_op_y*phi + (ny > 0) .* downwind_op_y*phi;

    advective_flux = velc*(nx .* dphi_dx + ny .* dphi_dy);


    phi = phi - dt*(advective_flux);

    rho = 1-smear_out_heaviside(phi,dx);


    % u = rho.*u;
    % v = rho.*v;

end



%%
vid = VideoWriter('UV_concentration_video9.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);

for i = 1:1000
    if mod(i,10) == 0
        clf

        subplot(2,1,1)
        z = reshape(Rho(:,i),grid_size,grid_size);
        contour(z, 'r', 'LineWidth', 5); colorbar; axis equal tight;
        hold on
        z = reshape(U(:,i),grid_size,grid_size);
        imagesc(z,'AlphaData', z ~= 0);colormap(parula); colorbar; axis equal tight;
        title('Concentration of U');


        subplot(2,1,2)

        z = reshape(Rho(:,i),grid_size,grid_size);
        contour(z, 'r', 'LineWidth', 5); colorbar; axis equal tight;
        hold on
        z = reshape(V(:,i),grid_size,grid_size);
        imagesc(z,'AlphaData', z ~= 0);colormap(parula); colorbar; axis equal tight;
        title('Concentration of U');

        frame = getframe(gcf);
        writeVideo(vid, frame);
        
    end
end 

close(vid);

%%

vid = VideoWriter('boundary2.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);

for i = 1:1000
    if mod(i,10) == 0
        clf

        z = reshape(Rho(:,i),grid_size,grid_size);
        imagesc(z(1:100,1:200));
        % contour(z, [0.5 0.5],'r', 'LineWidth', 2); colorbar; axis equal tight;

        frame = getframe(gcf);
        writeVideo(vid, frame);
        
    end
end 

close(vid);