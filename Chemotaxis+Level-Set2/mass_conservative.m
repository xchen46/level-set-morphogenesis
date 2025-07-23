%% Minimal Partial Differential Equation Model for Level-Set + Chemotaxis

% Version 3.0
% NOTE: Needs to fix the boundaer

%% Clear workspace and figures
clc
clear
close all

%% Initialize Grid Parameters
tic
membrane_size = 50;                  % Physical size of the domain
grid_size     = 258;                 % Number of grid points along one axis
num           = grid_size^2;        % Total number of grid entries
dx            = membrane_size / (grid_size - 1);  % Grid spacing
x             = linspace(0, membrane_size, grid_size); % Grid coordinates
[dom_x, dom_y] = meshgrid(x, x);    % Create 2D meshgrid

%% Simulation Parameters
dt   = 0.0005;           % Time step for evolution
dtau = dx * 0.005;       % Time step for reinitialization

Du   = 0.5;              % Diffusion coefficient of U
Dv   = 0.5;              % Diffusion coefficient of V
velc = 1;                % Advection velocity coefficient

a = 0.5;                 % Coupling term coefficient
b = 0.2;                 % Reaction term coefficient
c = 2;                   % Chemotaxis sensitivity coefficient

%% Initialize Fields
U = zeros(num, 1000);    % Allocate U field over time
V = zeros(num, 1000);    % Allocate V field over time


%% Initialize Fields: U, V, and Buffer Mask

% Create a logical mask to exclude buffer boundaries (outermost grid cells)
inside_mask = true(grid_size, grid_size);
inside_mask([1 end], :) = false;       % Exclude top and bottom rows
inside_mask(:, [1 end]) = false;       % Exclude left and right columns
inside_mask = inside_mask(:);          % Flatten mask to column vector

% Set random seed for reproducibility
rng(1)

% Initialize U with sinusoidal pattern in X-direction (λ ≈ 25)
f = sin(2 * pi * dom_x / 25);
u = 0.2 + 0.05 * f;                   % Small oscillation around 0.2
u = u(:);                             % Flatten to column vector
u(~inside_mask) = 0;                  % Zero out buffer boundary
U(:,1) = u;                           % Assign as initial condition

% Initialize V with sinusoidal pattern in Y-direction (λ ≈ 25)
f = sin(2 * pi * dom_y / 25);
v = 0.2 + 0.05 * f;
v = v(:);
v(~inside_mask) = 0;
V(:,1) = v;

%% Operator Construction
% Build all spatial differential operators and patch them to enforce Neumann boundary conditions.

% === Neumann boundary patches ===
% Different patching modes for testing boundary behavior
[Ops.neumann_patch_x, Ops.neumann_patch_y]   = apply_Neumann(sparse(num,num), grid_size, grid_size, 'one_sided');
[Ops.symmetric_patch_x, Ops.symmetric_patch_y] = apply_Neumann(sparse(num,num), grid_size, grid_size, 'symmetric');
[Ops.weird_x, Ops.weird_y]                   = apply_Neumann(sparse(num,num), grid_size, grid_size, 'weird');
% [Ops.weirder_x, Ops.weirder_y]             = apply_Neumann(sparse(num,num), grid_size, grid_size, 'weirder'); % unused

% Combine symmetric patch terms
Ops.symmetric_patch = Ops.symmetric_patch_x + Ops.symmetric_patch_y;

% === Identity matrix ===
Ops.speye = speye(num);

% === Laplacian (internal) operators ===
[Ops.inlap_x, Ops.inlap_y] = interial_laplace(grid_size, grid_size);


% === Central difference divergence operators ===
[Dx, Dy]             = interial_divergence_CDM(grid_size, grid_size);
[Ops.Dx, Ops.Dy]     = exposed_divergence(grid_size, grid_size);  % edge-aware version

% Add boundary patching (layer 2) to central difference derivatives
Ops.Dx2 = Dx / (2*dx) + Ops.neumann_patch_x / dx;
Ops.Dy2 = Dy / (2*dx) + Ops.neumann_patch_y / dx;

% Exposed divergence scaled (used in level-set normal estimation etc.)
Ops.Dx3 = Ops.Dx / (2*dx);  % Note: sign-sensitive
Ops.Dy3 = Ops.Dy / (2*dx);

% === Forward and Backward finite differences ===
[FWx, BWx] = interial_divergence_FDM(grid_size, grid_size, 1);  % x-direction
[FWy, BWy] = interial_divergence_FDM(grid_size, grid_size, 2);  % y-direction


% === Upwind and Downwind Operators (internal) ===
[Ops.upwind_op_x, Ops.downwind_op_x] = interial_updownwind(grid_size, grid_size, 1);
[Ops.upwind_op_y, Ops.downwind_op_y] = interial_updownwind(grid_size, grid_size, 2);

% === Upwind and Downwind Operators (exposed edges) ===
[Ops.ex_upwind_op_x, Ops.ex_downwind_op_x] = exposed_updownwind(grid_size, grid_size, 1);
[Ops.ex_upwind_op_y, Ops.ex_downwind_op_y] = exposed_updownwind(grid_size, grid_size, 2);

% === Scale up/downwind operators by dx ===
Ops.upwind_op_x    = Ops.upwind_op_x    / dx;
Ops.upwind_op_y    = Ops.upwind_op_y    / dx;
Ops.downwind_op_x  = Ops.downwind_op_x  / dx;
Ops.downwind_op_y  = Ops.downwind_op_y  / dx;

% === Optional/Experimental logic (commented for reference) ===
% These versions patch the exposed operators — only use one patching method at a time:
% Ops.ex_upwind_op_x = (Ops.ex_upwind_op_x - Ops.weird_x) / dx;
% Ops.ex_upwind_op_y = (Ops.ex_upwind_op_y - Ops.weird_y) / dx;
% Ops.ex_downwind_op_x = (Ops.ex_downwind_op_x + Ops.weird_x) / dx;
% Ops.ex_downwind_op_y = (Ops.ex_downwind_op_y + Ops.weird_y) / dx;

% === Notes ===
% Keep internal up/downwind operators as they generate internal flux.
% Interesting behavior: patching only one of upwind/downwind works, but not both.



%% Initialize Level-Set Function φ and Density Field ρ


% === Level-set initialization ===
phi = dom_y - 25;           % Initialize φ as a vertical interface at y = 25
% Alternative (commented): circular interface
% phi = sqrt((dom_x - 25).^2 + (dom_y - 25).^2) - 10;

phi = phi(:);               % Flatten φ to column vector
phi(~inside_mask) = 0;      % Zero out values outside the domain mask

% === Smoothed Heaviside and Delta functions ===
Rho        = zeros(num, 1000);            % Initialize ρ matrix
Rho(:, 1)  = 1 - smear_out_heaviside(phi, dx);  % Initial density field based on φ
del        = smear_out_delta(phi, dx);          % Smoothed δ-function of φ

% === Mask initial u, v fields by ρ ===
rho = Rho(:,1);             % Extract initial density

u = rho .* u;               % Apply density mask to u
u = u ./ sum(u);            % Normalize u to preserve total mass
v = rho .* v;               % Apply density mask to v


%% Time Evolution Loop
tend = 500;

for i = 2:tend

    % --- 1. Build Diffusion Matrices (Implicit) for u and v ---
    Dr = spdiags(rho, 0, num, num);  % Diagonal matrix of ρ
    % Dr = Ops.speye;                % Optional: disable density scaling for testing

    % Construct implicit diffusion operators with density-weighted Laplacian
    lapU = Ops.speye - Du * dt * (BWx * Dr * FWx + BWy * Dr * FWy + Ops.symmetric_patch) / dx^2;
    lapV = Ops.speye - Dv * dt * (BWx * Dr * FWx + BWy * Dr * FWy + Ops.symmetric_patch) / dx^2;

    % --- 2. Compute Chemotactic Flux ---
    % Determine flux directions for v
    flux_sign_x = sign(Ops.inlap_x * v);
    flux_sign_y = sign(Ops.inlap_y * v);

    % Directional gradients of v using upwind/downwind
    grad_v_x = (flux_sign_x >= 0) .* (Ops.upwind_op_x * v) + ...
               (flux_sign_x <  0) .* (Ops.downwind_op_x * v);

    grad_v_y = (flux_sign_y >= 0) .* (Ops.upwind_op_y * v) + ...
               (flux_sign_y <  0) .* (Ops.downwind_op_y * v);

    % Compute chemotactic flux J = ρ * sensitivity * ∇v
    sensitivity = c * u ./ (1 + u.^2 + 1e-3);  % Chemotactic sensitivity
    Jx = rho .* (sensitivity .* grad_v_x);
    Jy = rho .* (sensitivity .* grad_v_y);

    % Compute divergence of chemotactic flux
    div_Jv = Ops.Dx2 * Jx + Ops.Dy2 * Jy;

    % --- 3. Compute Reaction Term for v ---
    fv = rho .* (a * u - b * v);

    % --- 4. Time Step u and v using Implicit Schemes ---
    % Solve for u and v implicitly (can use GMRES if needed)
    u = lapU \ (u - dt * div_Jv);   % Advection–diffusion
    v = lapV \ (v + dt * fv);      % Reaction–diffusion

    % --- 5. Save Snapshot Every 20 Steps ---
    if mod(i, 20) == 0
        U(:, i/20)   = u;
        V(:, i/20)   = v;
        Rho(:, i/20) = phi;

        % (Optional) Reinitialization of φ to preserve signed-distance property
        % phi0 = phi;
        % sign_phi0 = sign(phi0);
        % for k = 1:10
        %     grad_x = (sign_phi0 >= 0) .* (Ops.ex_upwind_op_x * phi) + ...
        %              (sign_phi0 <  0) .* (Ops.ex_downwind_op_x * phi);
        %     grad_y = (sign_phi0 >= 0) .* (Ops.ex_upwind_op_y * phi) + ...
        %              (sign_phi0 <  0) .* (Ops.ex_downwind_op_y * phi);
        %     norm_grad_phi = sqrt(grad_x.^2 + grad_y.^2);
        %     phi = phi - dtau * sign(sign_phi0) .* (norm_grad_phi - 1);
        % end
    end

    % --- 6. Evolve Level-Set φ via Advection ---
    % Compute normal direction of φ
    [nx, ny] = phi_normal(phi, Ops.Dx3, Ops.Dy3);

    % Select upwind derivatives based on normal direction
    sign_nx = sign(nx); sign_ny = sign(ny);
    dphi_dx = (sign_nx >= 0) .* (Ops.ex_upwind_op_x * phi) + ...
              (sign_nx <  0) .* (Ops.ex_downwind_op_x * phi);
    dphi_dy = (sign_ny >= 0) .* (Ops.ex_upwind_op_y * phi) + ...
              (sign_ny <  0) .* (Ops.ex_downwind_op_y * phi);

    % Advection update for φ
    advective_flux = velc * (nx .* dphi_dx + ny .* dphi_dy);
    phi = phi - dt * advective_flux;

    % --- 7. Update Density Field ρ Based on φ ---
    rho = 1 - smear_out_heaviside(phi, dx);
end

toc

%% Visualization and Video Export (Optional)
% Create a video showing the time evolution of U and V

% % === Set up video writer ===
vid = VideoWriter('moving_bd.mp4', 'MPEG-4');  % Create video file
vid.FrameRate = 10;                            % Set frame rate
open(vid);                                     % Open file for writing

% % === Loop through saved snapshots ===
for i = 1:(tend / 20)  % Match snapshot saving frequency in main loop
    if mod(i, 1) == 0  % Plot every frame (can skip frames if needed)
        clf

        % --- Plot U concentration ---
        subplot(2,1,1)
        z = reshape(U(:,i), grid_size, grid_size);  % Reshape U to 2D
        hU = imagesc(z);                            % Plot as image
        set(hU, 'AlphaData', z ~= 0);               % Mask zeros (optional)
        colormap(viridis); colorbar;
        axis equal tight;
        title('Concentration of U');

        % --- Plot V concentration ---
        subplot(2,1,2)
        z = reshape(V(:,i), grid_size, grid_size);  % Reshape V to 2D
        hV = imagesc(z);                            % Plot as image
        set(hV, 'AlphaData', z ~= 0);               % Mask zeros (optional)
        colormap(viridis); colorbar;
        axis equal tight;
        title('Concentration of V');

        % --- Save current frame to video ---
        frame = getframe(gcf);        % Capture current figure
        writeVideo(vid, frame);       % Write frame to video
    end
end

% % === Close video writer ===
close(vid);

%% Video Output: Boundary Evolution of φ / ρ (fixing needed)
% This script generates a video showing the evolution of Rho (based on φ) over time.

vid = VideoWriter('boundary_evolve.mp4', 'MPEG-4');  % Output file name and format
vid.FrameRate = 10;                            % Playback speed (frames per second)
open(vid);                                     % Open video file for writing

for i = 1:(tend / 20)  % Iterate through saved snapshots (saved every 20 steps)
    if mod(i, 20) == 0  % Optional frame-skipping logic; currently redundant since i steps in 20s
        clf  % Clear current figure

        % --- Extract and reshape Rho for current frame ---
        z = reshape(Rho(:, i), grid_size, grid_size);  % Convert vector to 2D field

        % --- Display density map ---
        imagesc(z);                % Visualize Rho using color intensity
        colorbar;
        colormap(viridis);         % Optional: apply colormap
        axis equal tight;          % Fix aspect ratio and remove excess padding
        title(sprintf('Rho at Step %d', i * 20));

        % Optional: Add φ = 0.5 contour line
        % contour(z, [0.5 0.5], 'r', 'LineWidth', 2);

        % --- Save frame to video ---
        frame = getframe(gcf);
        writeVideo(vid, frame);
    end
end

close(vid);  % Finalize video file


%% Boundary Behavior Visualization (Optional)
% This script generates a video showing differences between inner and outer boundary rows/columns of V

% % === Set up video writer ===
vid = VideoWriter('neuman_bd.mp4', 'MPEG-4');  % Output video file
vid.FrameRate = 10;                               % Frames per second
open(vid);                                        % Open file for writing

% % === Loop through snapshots ===
for i = 1:(tend / 20)  % Match snapshot frequency
    if mod(i, 1) == 0
        clf

        % --- Reshape V back to 2D for this time step ---
        z = reshape(V(:, i), grid_size, grid_size);

        % --- Plot boundary differences to examine behavior ---
        % Differences in y-direction (top and bottom)
        plot(z(3, 2:end-1) - z(2, 2:end-1), 'r');     % Near top boundary
        hold on
        plot(z(end-2, 2:end-1) - z(end-1, 2:end-1), 'b');  % Near bottom boundary

        % Differences in x-direction (left and right)
        plot(z(2:end-1, 3) - z(2:end-1, 2), 'g');     % Near left boundary
        plot(z(:, end-2) - z(:, end-1), 'm'); % Near right boundary

        title('Boundary Gradient Check: V');
        % legend({'Top', 'Bottom', 'Left', 'Right'}, 'Location', 'best');

        % --- Save frame to video ---
        frame = getframe(gcf);
        writeVideo(vid, frame);
    end
end

% % === Close video writer ===
close(vid);