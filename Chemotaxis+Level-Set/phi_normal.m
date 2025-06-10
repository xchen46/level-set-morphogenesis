
function [grad_phi_x,grad_phi_y] = phi_normal(phi,Dx,Dy)

    grad_phi_x = Dx*phi;
    grad_phi_y = Dy*phi;
    norm = sqrt(grad_phi_x.^2+grad_phi_y.^2 + 1e-10);

    grad_phi_x = grad_phi_x./norm;
    grad_phi_y = grad_phi_y./norm;

end
