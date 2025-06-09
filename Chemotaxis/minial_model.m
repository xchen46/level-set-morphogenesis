%% Minimal Partial Differenital Equation Model for Chemotaxis%% 
%%
clc
clear 
close all
%%
membrane_size = 5;
grid_size = 200;
dx = membrane_size/(grid_size - 1); 
x = linspace(0,membrane_size,grid_size);
[dom_x,dom_y]  = meshgrid(x,x);
%%

u = bolb_rand(x,membrane_size);
[v,~]  = meshgrid(linspace(0,0.1,grid_size),linspace(0,0.1,grid_size));
%%

% cx1 = 2; cy1 = 2.5;   % Left blob
% cx2 = 3; cy2 = 2.5;   % Right blob
% sigma_u = 0.3;
% 
% % Create two Gaussians
% u1 = exp(-((dom_x - cx1).^2 + (dom_y - cy1).^2) / (2*sigma_u^2));
% u2 = exp(-((dom_x - cx2).^2 + (dom_y - cy2).^2) / (2*sigma_u^2));
% u = u1 + u2;
% 
% % Initialize v as a smoothed version of u
% sigma_v = 1.5;
% kernel_size = 7;
% [xk, yk] = meshgrid(-kernel_size:kernel_size, -kernel_size:kernel_size);
% G = exp(-(xk.^2 + yk.^2) / (2 * sigma_v^2));
% G = G / sum(G(:));
% v = conv2(u, G, 'same')/10;
%%
% v = exp(-((dom_x - 1.5).^2 + (dom_y - 2.5).^2)/0.1) + ...
%     0.5 * exp(-((dom_x - 3).^2 + (dom_y - 2.5).^2)/0.2);
% v = zeros(size(dom_x));
dt = 0.0001;
D = 0.2;
kappa = 1;
%% nonflux boundary initialize
u(1,:) = u(3,:);
u(2,:) = u(3,:);
u(grid_size,:) = u(grid_size-2,:);
u(grid_size-1,:) = u(grid_size-2,:);

u(:,1) = u(:,3);
u(:,2) = u(:,3);
u(:,grid_size) = u(:,grid_size-2);
u(:,grid_size-1) = u(:,grid_size-2);

v(1,:) = v(3,:);
v(2,:) = v(3,:);
v(grid_size,:) = v(grid_size-1,:);
v(grid_size-1,:) = v(grid_size-1,:);

v(:,1) = v(:,3);
v(:,2) = v(:,3);
v(:,grid_size) = v(:,grid_size-2);
v(:,grid_size-1) = v(:,grid_size-2);


%%
U = zeros(grid_size,grid_size,5000);
V = zeros(grid_size,grid_size,5000);
U(:,:,1) = u;
V(:,:,1) = v;

u = zeros(grid_size);
v = zeros(grid_size);
flux_x = zeros(grid_size);
flux_y = zeros(grid_size);
u_next = u;
v_next = v;
%%
N = grid_size - 2;


for i = 2:10000
    u = U(:,:,i-1);
    v = V(:,:,i-1);
    
    for j = 1:N^2
        a = floor((j-1)/N)+2;
        b = mod(j-1,N)+2; 

        up_x = (kappa * (v(a+1,b) - v(a-1,b))) >= 0;
        up_y = (kappa * (v(a,b+1) - v(a,b-1))) >= 0;

        flux_x(a,b) = -D*upwind(u,dx,a,b,0,0)+kappa*u(a,b)*upwind(v,dx,a,b,up_x,0);
        flux_y(a,b) = -D*upwind(u,dx,a,b,0,1)+kappa*u(a,b)*upwind(v,dx,a,b,up_y,1);

        
    end

    for k = 1:N^2
        a = floor((k-1)/N)+2;
        b = mod(k-1,N)+2; 
        
        u_next(a,b) = -(flux_x(a+1,b)-flux_x(a-1,b))/(2*dx)-...
            (flux_y(a,b+1)-flux_y(a,b-1))/(2*dx);
        v_next(a,b) = (CDM(v,dx,a,b,0)+ CDM(v,dx,a,b,1))+...
                 u(a,b)-v(a,b);
    end
    
    U(:,:,i) = u+dt*u_next;
    V(:,:,i) = v+dt*v_next;    
    
    U(1,:,i) = U(3,:,i);
    U(grid_size,:,i) = U(grid_size-2,:,i);
    U(:,1,i) = U(:,3,i);
    U(:,grid_size,i) = U(:,grid_size-2,i);

    V(1,:,i) = V(3,:,i);
    V(grid_size,:,i) = V(grid_size-2,:,i);
    V(:,1,i) = V(:,3,i);
    V(:,grid_size,i) = V(:,grid_size-2,i);
end



%%
vid = VideoWriter('UV_concentration_video5.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);

for i = 1:10000
    if mod(i,40) == 0
        clf
        subplot(2,1,1)
        imagesc(U(:,:,i)); colorbar; axis equal tight;
        title('Concentration of U');

        subplot(2,1,2)
        imagesc(V(:,:,i)); colorbar; axis equal tight;
        title('Concentration of V');

        frame = getframe(gcf);
        writeVideo(vid, frame);
        
    end
end 

close(vid);

%%
vid = VideoWriter('UV_gradient_video.mp4', 'MPEG-4');
vid.FrameRate = 10;  % Adjust frame rate as needed
open(vid);

for i = 1:10000
    if mod(i,40) == 0
        clf
        [dvx, dvy] = gradient(U(:,:,i), 600*dx);
        quiver(dom_x, dom_y, dvx, dvy,6);
        axis equal tight;
        title('Concentration Gradient of U');

        
        frame = getframe(gcf);
        writeVideo(vid, frame);
    end
end 

close(vid);
