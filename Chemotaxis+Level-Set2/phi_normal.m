function [grad_phi_x, grad_phi_y] = phi_normal(phi, Dx, Dy)
%PHI_NORMAL Compute the normalized gradient (unit normal) of phi
%
%   Inputs:
%       phi - column vector representing the scalar field (flattened 2D grid)
%       Dx  - sparse matrix for discrete ∂/∂x operator
%       Dy  - sparse matrix for discrete ∂/∂y operator
%
%   Outputs:
%       grad_phi_x - x-component of the normalized gradient (unit vector)
%       grad_phi_y - y-component of the normalized gradient (unit vector)

    % Compute gradient components using finite difference operators
    grad_phi_x = Dx * phi;
    grad_phi_y = Dy * phi;

    % Compute the magnitude of the gradient (add small ε to avoid divide-by-zero)
    norm = sqrt(grad_phi_x.^2 + grad_phi_y.^2 + 1e-10);

    % Normalize the gradient to obtain unit normal vectors
    grad_phi_x = grad_phi_x ./ norm;
    grad_phi_y = grad_phi_y ./ norm;
end

