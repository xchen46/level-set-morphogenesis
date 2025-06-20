%% Initialize Parameters
membrane_size = 50; % Membrane Size
grid_size = 250; % Number of Grid_Point
num = grid_size^2; % Total Number of entries
dx = membrane_size/(grid_size - 1);  % Grid step
x = linspace(0,membrane_size,grid_size); 
[dom_x,dom_y]  = meshgrid(x,x); % Create grid
%% Operator Consturction
% Need to make sure the boundary conditions is doing its job

Ops.neumann_patch = apply_Neumann(sparse(num,num),grid_size,grid_size,'one_sided');
Ops.symmetric_patch = apply_Neumann(sparse(num,num),grid_size,grid_size,'symmetric');

[Ops.inlap_x, Ops.inlap_y] = interial_laplace(grid_size,grid_size);
Ops.lapc = (Ops.inlap_x+Ops.inlap_y)./(2*dx^2)+Ops.symmetric_patch./(2*dx^2);

[Dx, Dy] = interial_divergence(grid_size,grid_size,dx,'unbound');
Ops.Dx = Dx + Ops.neumann_patch./dx;
Ops.Dy = Dy + Ops.neumann_patch./dx;

[Ops.upwind_op_x, Ops.downwind_op_x] = interial_updownwind(grid_size,grid_size,1);
[Ops.upwind_op_y, Ops.downwind_op_y] = interial_updownwind(grid_size,grid_size,2);

Ops.upwind_op_x = (Ops.upwind_op_x+Ops.neumann_patch)./dx;
Ops.upwind_op_y = (Ops.upwind_op_y+Ops.neumann_patch)./dx;
Ops.downwind_op_x = (Ops.downwind_op_x+Ops.neumann_patch)./dx;
Ops.downwind_op_y = (Ops.downwind_op_y+Ops.neumann_patch)./dx;


%% Simulation Parameters

dt = 0.005;
D = 1;
% kappa = 0.3;
velc = 0.1;
a = 0.1;
c = 4;

U = zeros(num,1000);
V = zeros(num,1000);

Rho = zeros(num,1000);
Rho(:,1) = 1-smear_out_heaviside(phi,dx);

u = Rho(:,1).*u;
v = Rho(:,1).*v;

