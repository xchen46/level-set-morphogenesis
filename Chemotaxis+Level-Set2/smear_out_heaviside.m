function res = smear_out_heaviside(phi, dx)
% SMEAR_OUT_HEAVISIDE
% Computes a smoothed Heaviside function based on the signed distance field φ.
% Useful in level-set methods for smoothly approximating interfaces.
%
% Inputs:
%   phi - signed distance function (vector or matrix)
%   dx  - grid spacing (used to define smoothing width)
%
% Output:
%   res - smoothed Heaviside values in [0, 1]

    epsilon = 1.5 * dx;  % Width of smoothing band

    res = zeros(size(phi));  % Preallocate output

    % --- Region: φ < -ε (fully inside the interface)
    idx = phi < -epsilon;
    res(idx) = 0;

    % --- Region: |φ| ≤ ε (transition zone)
    idx = (-epsilon < phi) & (phi < epsilon);
    res(idx) = 1/2 + phi(idx) / (2 * epsilon) + ...
               1/(2*pi) * sin(pi * phi(idx) / epsilon);

    % --- Region: φ > ε (fully outside the interface)
    idx = phi > epsilon;
    res(idx) = 1;
end
