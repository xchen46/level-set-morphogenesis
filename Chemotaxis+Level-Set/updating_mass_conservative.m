%% Minimal Partial Differenital Equation Model for Level-Set + Chemotaxis%% 

% update 2.0
% changed to implicit diffusion but the boundary is off now
%% clear all the workspace, close all the windows, and clear all command
clc
clear 
close all
%% Initialize Parameters
tic
membrane_size = 50; % Membrane Size
grid_size = 256; % Number of Grid_Point
num = grid_size^2; % Total Number of entries
dx = membrane_size/(grid_size - 1);  % Grid step
x = linspace(0,membrane_size,grid_size); 
[dom_x,dom_y]  = meshgrid(x,x); % Create grid

%% Simulation Parameters

dt = 0.0005;
dtau = dx*0.005; % reinitialization step
% D = 2;
% % kappa = 0.3;
% velc = 0.1;
% a = 0.5;
% c = 0.5;

Du = 0.5;
Dv = 0.5;
velc = 1;
a = 0.5;
b = 0.2;
c = 1;

U = zeros(num,1000);
V = zeros(num,1000);



%% Initialize U,V,phi

% [u,~] = initalize_two_blob(dom_x,dom_y,10,30,15,15); % initialize u and v
% [v,~] = initalize_two_blob(dom_x,dom_y,13,33,15,15); % initialize u and v


% [xg, yg] = meshgrid(-5:5, -5:5);
% gauss_kernel = exp(-(xg.^2 + yg.^2)/(2*2^2));
% gauss_kernel = gauss_kernel / sum(gauss_kernel(:));
% Smooth the initial field

rng(1)

f = sin(2*pi*dom_x/25);  % wavelength ~25
u = 0.2 + 0.05 * f;
u = u(:);


% u = 0.2+0.1*randn(grid_size,grid_size);

% u = max(u, 0);
% u = u(:);

% v = 0.2+0.1*randn(grid_size,grid_size);

f = sin(2*pi*dom_y/25);  % wavelength ~25
v = 0.2 + 0.05 * f;
v = v(:);



% v = max(v, 0);
% v = v(:);


%%

interior_mask = true(grid_size, grid_size);
interior_mask([1 end], :) = false;
interior_mask(:, [1 end]) = false;
interior_mask = interior_mask(:);


%% Operator Consturction
% Need to make sure the boundary conditions is doing its job

[Ops.neumann_patch_y,Ops.neumann_patch_x] = apply_Neumann(sparse(num,num),grid_size,grid_size,'one_sided');
[Ops.symmetric_patch,~] = apply_Neumann(sparse(num,num),grid_size,grid_size,'symmetric');

[Ops.inlap_x, Ops.inlap_y] = interial_laplace(grid_size,grid_size);
Ops.lapc = (Ops.inlap_x+Ops.inlap_y)./(dx^2)+2*Ops.symmetric_patch./(dx^2);
% Ops.symmetric_patch = 2*Ops.symmetric_patch./(dx^2);

Ops.speye = speye(num);
% Ops.lapcU = Ops.speye - dt*Ops.lapc;
% Ops.lapcV = Ops.speye - dt*D*Ops.lapc;

[Dx, Dy] = interial_divergence(grid_size,grid_size,dx,'unbound');
Ops.Dx = Dx;
Ops.Dy = Dy;
Ops.Dx2 = Dx + (Ops.neumann_patch_x)./(dx);
Ops.Dy2 = Dy + (Ops.neumann_patch_y)./(dx);

[Dx3, Dy3] = full_divergence(grid_size,grid_size,dx,'unbound');

[Ops.upwind_op_x, Ops.downwind_op_x] = interial_updownwind(grid_size,grid_size,1);
[Ops.upwind_op_y, Ops.downwind_op_y] = interial_updownwind(grid_size,grid_size,2);

% [Ops.upwind_op_x2, Ops.downwind_op_x2] = exposed_updownwind(grid_size,grid_size,1);
% [Ops.upwind_op_y2, Ops.downwind_op_y2] = exposed_updownwind(grid_size,grid_size,2);


Ops.upwind_op_x = (Ops.upwind_op_x)./dx;
Ops.upwind_op_y = (Ops.upwind_op_y)./dx;
Ops.downwind_op_x = (Ops.downwind_op_x)./dx;
Ops.downwind_op_y = (Ops.downwind_op_y)./dx;

% Ops.upwind_op_x2 = (Ops.upwind_op_x2)./dx;
% Ops.upwind_op_y2 = (Ops.upwind_op_y2)./dx;
% Ops.downwind_op_x2 = (Ops.downwind_op_x2)./dx;
% Ops.downwind_op_y2 = (Ops.downwind_op_y2)./dx;

% keep the upwind+downwind internal since they give us flux which should be 0 
% it is interesting if you patch one operator or the other it works, but
% not both

%%
% phi_mask = true(grid_size,grid_size);  % assume phi is reshaped to 2D
% phi_mask(1, :) = false;      % top
% phi_mask(end, :) = false;    % bottom
% phi_mask(:, 1) = false;      % left
% phi_mask(:, end) = false;    % right
% phi_mask = phi_mask(:);

phi = dom_y-25;
% phi = sqrt((dom_x - 25).^2 + (dom_y - 25).^2)-10;
phi = phi(:);
Rho = zeros(num,1000);
Rho(:,1) = 1-smear_out_heaviside(phi,dx);
del = smear_out_delta(phi,dx);

Rho(~interior_mask,1) = 0;
rho = Rho(:,1);


u = rho.*u;
v = rho.*v;
%%
tend = 10000;
for i = 2:tend
    
    Dr  = spdiags(rho,0,num,num);

    lapU  = Ops.speye-Du*dt*(Ops.Dx*Dr*Ops.Dx + Ops.Dy*Dr*Ops.Dy + Dr*2*Ops.symmetric_patch./(dx^2));
    lapV  = Ops.speye-Dv*dt*(Ops.Dx*Dr*Ops.Dx + Ops.Dy*Dr*Ops.Dy + Dr*2*Ops.symmetric_patch./(dx^2));

    % fprintf('Condition estimate of lapU: %.2e\n', condest(lapU));

    % % Simulating what is inside the boundary
    flux_sign_x = sign(Ops.inlap_x*v);
    flux_sign_y = sign(Ops.inlap_y*v);
    grad_v_x = (flux_sign_x >= 0) .* Ops.upwind_op_x*v + (flux_sign_x < 0) .* Ops.downwind_op_x*v;
    grad_v_y = (flux_sign_y >= 0) .* Ops.upwind_op_y*v + (flux_sign_y < 0) .* Ops.downwind_op_y*v;


    Jx = rho.*(c*u./(1+u.^2+ 1e-3).*grad_v_x);
    Jy = rho.*(c*u./(1+u.^2+ 1e-3).*grad_v_y);

    div_Jv = Ops.Dx2 *Jx+ Ops.Dy2 *Jy;

    fv = rho.*(a*u-b*v);

    [uu,~] = gmres(lapU,u, [], 1e-10, 200);
    % now i think there is a porblem in gmres
    % u = lapU \ (u-dt*div_Jv);
    % u = lapU \ (u);
    u = dt*(-div_Jv) + uu;

    [vv,~] = gmres(lapV,v, [], 1e-10, 200);
    % v = lapV \ (v+dt*fv);
    % vv = lapV \ v;
    v = dt*(fv) + vv;



    if mod(i,20) == 0
        U(:,i/20) = u;
        V(:,i/20) = v;
        Rho(:,i/20) = phi;

        
        phi0 = phi;  % Freeze the pre-evolved φ for reinit
        sign_phi0 = sign(phi0);

        % for k = 1:15  % evolve in τ
        %     grad_x = (sign_phi0 >= 0) .* (Ops.upwind_op_x * phi) + ...
        %         (sign_phi0 <  0) .* (Ops.downwind_op_x * phi);
        % 
        %     grad_y = (sign_phi0 >= 0) .* (Ops.upwind_op_y * phi) + ...
        %         (sign_phi0 <  0) .* (Ops.downwind_op_y * phi);
        % 
        %     norm_grad_phi = sqrt(grad_x.^2 + grad_y.^2);
        % 
        %     phi = phi - dtau * sign(sign_phi0) .* (norm_grad_phi - 1);
        % end
    end

    % Simulating rho
    % [nx, ny] = phi_normal(phi,Dx,Dy);
    % sign_nx = sign(nx); sign_ny = sign(ny);
    % dphi_dx = (sign_nx >= 0) .* Ops.upwind_op_x2*phi + (sign_nx < 0) .* Ops.downwind_op_x2*phi;
    % dphi_dy = (sign_ny >= 0) .* Ops.upwind_op_y2*phi + (sign_ny < 0) .* Ops.downwind_op_y2*phi;
    % 
    % advective_flux = velc*(nx .* dphi_dx + ny .* dphi_dy);
    % 
    % advective_flux(~interior_mask) = 0;
    % phi = phi - dt*advective_flux;
    % 
    % rho = 1-smear_out_heaviside(phi,dx);
    % rho(~interior_mask) = 0;
   
end

toc

%%
vid = VideoWriter('moving_bd6.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);
% clims = [0 1];
for i = 1:tend/20
    if mod(i,1) == 0
        clf

        subplot(2,1,1)
        % z = reshape(Rho(:,i),grid_size,grid_size);
        % contour(z, 'r', 'LineWidth', 5); axis equal tight;
        % hold on
        z = reshape(U(:,i),grid_size,grid_size);
        k = imagesc(z); 
        set(k, 'AlphaData', z ~= 0);
        colormap(viridis); colorbar; axis equal tight;
        title('Concentration of U');


        subplot(2,1,2)

        % z = reshape(Rho(:,i),grid_size,grid_size);
        % contour(z, 'r', 'LineWidth', 5); axis equal tight;
        % hold on
        z = reshape(V(:,i),grid_size,grid_size);
        k = imagesc(z); 
        set(k, 'AlphaData', z ~= 0);
        colormap(viridis); colorbar; axis equal tight;
        title('Concentration of V');

        frame = getframe(gcf);
        writeVideo(vid, frame);

    end
end 

close(vid);
 
%%
% 
vid = VideoWriter('comparsion6bd.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);
clims = [0 1];
for i = 1:tend/20
    if mod(i,20) == 0
        clf

        z = reshape(V(:,i),grid_size,grid_size);

        % imagesc(z(:,end-1:end)), colorbar;
        plot(z(2,:)-z(1,:));
        hold on

        plot(z(end-1,:)-z(end,:));

        plot(z(:,2)-z(:,1));
        plot(z(:,end-1)-z(:,end));

        title('Concentration of U');

        frame = getframe(gcf);
        writeVideo(vid, frame);

    end
end 

close(vid);



%% 
vid = VideoWriter('boundary6.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);

for i = 1:tend/20
    if mod(i,1) == 0
        clf

        z = reshape(Rho(:,i),grid_size,grid_size);
        imagesc(z);
        % contour(z, [0.5 0.5],'r', 'LineWidth', 2); colorbar; axis equal tight;

        frame = getframe(gcf);
        writeVideo(vid, frame);

    end
end 

close(vid);
%%

% vid = VideoWriter('UV_edge2.mp4', 'MPEG-4');
% vid.FrameRate = 10;
% open(vid);
% for i = 1:1000
% clf
% z = U(251:500, i)-U(1:250, i);  % make sure U is [100 × 1000] or similar
% k = plot(z);
% colormap(viridis); colorbar;
% frame = getframe(gcf);
% writeVideo(vid, frame);
% end
% close(vid);