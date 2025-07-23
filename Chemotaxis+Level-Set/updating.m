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

dt = 0.005;
% D = 2;
% % kappa = 0.3;
% velc = 0.1;
% a = 0.5;
% c = 0.5;

D = 0.1;
velc = 0.1;
a = 0.5;
c = 2.3;

% U = zeros(num,1000);
% V = zeros(num,1000);
U = zeros(num,1);
V = zeros(num,1);

% Rho = zeros(num,1000);
% Rho(:,1) = 1-smear_out_heaviside(phi,dx);
% 
% u = Rho(:,1).*u;
% v = Rho(:,1).*v;

%% Initialize U,V,phi

% [u,v] = initalize_two_blob(dom_x,dom_y,4,6,1.5,1.5); % initialize u and v


[xg, yg] = meshgrid(-5:5, -5:5);
gauss_kernel = exp(-(xg.^2 + yg.^2)/(2*2^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel(:));
% Smooth the initial field

rng(1)
u = 1+randn(grid_size,grid_size)*0.01;
% u = conv2(u, gauss_kernel, 'same');


% u(1,:)     = 1;
% u(end,:)   = 1;
% u(:,1)     = 1;
% u(:,end)   = 1;

% u(1,:) = u(2,:);
% u(grid_size,:) = u(grid_size-1,:);
% u(:,1) = u(:,2);
% u(:,grid_size) = u(:,grid_size-1);

u = u(:);

v = 1/a+randn(grid_size,grid_size)*0.01;
% v = conv2(v, gauss_kernel, 'same');

% v(1,:) = v(2,:);
% v(grid_size,:) = v(grid_size-1,:);
% v(:,1) = v(:,2);
% v(:,grid_size) = v(:,grid_size-1);


% v(1,:)     = 10;
% v(end,:)   = 10;
% v(:,1)     = 10;
% v(:,end)   = 10;


v = v(:);

phi = dom_y-25;
phi = phi(:);

%% Operator Consturction
% Need to make sure the boundary conditions is doing its job

[Ops.neumann_patch_y,Ops.neumann_patch_x] = apply_Neumann(sparse(num,num),grid_size,grid_size,'one_sided');
[Ops.symmetric_patch,~] = apply_Neumann(sparse(num,num),grid_size,grid_size,'symmetric');

[Ops.inlap_x, Ops.inlap_y] = interial_laplace(grid_size,grid_size);
Ops.lapc = (Ops.inlap_x+Ops.inlap_y)./(dx^2)+2*Ops.symmetric_patch./(dx^2);

Ops.speye = speye(num);
Ops.lapcU = Ops.speye - dt*Ops.lapc;
Ops.lapcV = Ops.speye - dt*D*Ops.lapc;

[Dx, Dy] = interial_divergence(grid_size,grid_size,dx,'unbound');
Ops.Dx = (Dx + Ops.neumann_patch_x)./dx;
Ops.Dy = (Dy + Ops.neumann_patch_y)./dx;

[Ops.upwind_op_x, Ops.downwind_op_x] = interial_updownwind(grid_size,grid_size,1);
[Ops.upwind_op_y, Ops.downwind_op_y] = interial_updownwind(grid_size,grid_size,2);
% keep the upwind+downwind internal since they give us flux which should be 0 
% it is interesting if you patch one operator or the other it works, but
% not both



%%

for i = 2:40000
    
    % Simulating what is inside the boundary
    flux_sign_x = sign(Ops.inlap_x*v);
    flux_sign_y = sign(Ops.inlap_y*v);
    grad_v_x = (flux_sign_x >= 0) .* Ops.upwind_op_x*v + (flux_sign_x < 0) .* Ops.downwind_op_x*v;
    grad_v_y = (flux_sign_y >= 0) .* Ops.upwind_op_y*v + (flux_sign_y < 0) .* Ops.downwind_op_y*v;

    % Jx = kappa * (u.* grad_v_x);
    % Jy = kappa * (u.* grad_v_y);


    Jx = c*u./(1+u.^2+ 1e-3).*grad_v_x;
    Jy = c*u./(1+u.^2+ 1e-3).*grad_v_y;

    div_Jv = Ops.Dx *Jx+ Ops.Dy *Jy;

    fu = u.*(1-u);
    fv = u-a*v;

    % [uu,~] = gmres(Ops.lapcU,u, [], 1e-10, 200);
    % now i think there is a porblem in gmres
    uu = Ops.lapcU \ u;
    u = dt*(-div_Jv+fu) + uu;

    % [vv,~] = gmres(Ops.lapcV,v, [], 1e-10, 200);
    vv = Ops.lapcV \ v;
    v = dt*(fv) + vv;


    if mod(i,40000) == 0
        U(:,1) = u;
        V(:,1) = v;
        % Rho(:,i/20) = rho;
    end


    % Simulating rho
    % [nx, ny] = phi_normal(phi,Ops.Dx,Ops.Dy);
    % dphi_dx = (nx < 0) .* Ops.upwind_op_x*phi + (nx > 0) .* Ops.downwind_op_x*phi;
    % dphi_dy = (ny < 0) .* Ops.upwind_op_y*phi + (ny > 0) .* Ops.downwind_op_y*phi;
    % 
    % advective_flux = velc*(nx .* dphi_dx + ny .* dphi_dy);
    % 
    % 
    % phi = phi - dt*(advective_flux);
    % 
    % rho = 1-smear_out_heaviside(phi,dx);
    % 
    % 
    % u = rho.*u;
    % v = rho.*v;

end

toc
save('comparsion1_c23_2.mat','U','V');

%%
% vid = VideoWriter('comparsion2.mp4', 'MPEG-4');
% vid.FrameRate = 10;  % Adjust frame rate as needed
% open(vid);
% % clims = [0 1];
% for i = 1:1000
%     if mod(i,1) == 0
%         clf
% 
%         subplot(2,1,1)
%         % z = reshape(Rho(:,i),grid_size,grid_size);
%         % contour(z, 'r', 'LineWidth', 5); axis equal tight;
%         % hold on
%         z = reshape(U(:,i),grid_size,grid_size);
%         k = imagesc(z); 
%         set(k, 'AlphaData', z ~= 0);
%         colormap(viridis); colorbar; axis equal tight;
%         title('Concentration of U');
% 
% 
%         subplot(2,1,2)
% 
%         % z = reshape(Rho(:,i),grid_size,grid_size);
%         % contour(z, 'r', 'LineWidth', 5); axis equal tight;
%         % hold on
%         z = reshape(V(:,i),grid_size,grid_size);
%         k = imagesc(z); 
%         set(k, 'AlphaData', z ~= 0);
%         colormap(viridis); colorbar; axis equal tight;
%         title('Concentration of V');
% 
%         frame = getframe(gcf);
%         writeVideo(vid, frame);
% 
%     end
% end 
% 
% close(vid);
% 
% %%
% 
% vid = VideoWriter('comparsion2bd.mp4', 'MPEG-4');
% vid.FrameRate = 10;  % Adjust frame rate as needed
% open(vid);
% clims = [0 1];
% for i = 1:1000
%     if mod(i,1) == 0
%         clf
% 
%         z = reshape(V(:,i),grid_size,grid_size);
% 
%         % imagesc(z(:,end-1:end)), colorbar;
%         plot(z(2,:)-z(1,:));
%         hold on
% 
%         plot(z(end-1,:)-z(end,:));
% 
%         plot(z(:,2)-z(:,1));
%         plot(z(:,end-1)-z(:,end));
% 
%         title('Concentration of U');
% 
%         frame = getframe(gcf);
%         writeVideo(vid, frame);
% 
%     end
% end 
% 
% close(vid);
% 


%% 
% vid = VideoWriter('boundary11.mp4', 'MPEG-4');
% vid.FrameRate = 10;  % Adjust frame rate as needed
% open(vid);
% 
% for i = 1:1000
%     if mod(i,1) == 0
%         clf
% 
%         z = reshape(Rho(:,i),grid_size,grid_size);
%         imagesc(z);
%         % contour(z, [0.5 0.5],'r', 'LineWidth', 2); colorbar; axis equal tight;
% 
%         frame = getframe(gcf);
%         writeVideo(vid, frame);
% 
%     end
% end 
% 
% close(vid);
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