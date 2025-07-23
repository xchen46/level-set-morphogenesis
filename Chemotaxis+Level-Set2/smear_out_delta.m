function res = smear_out_delta(phi, dx)
% SMEAR_OUT_DELTA
% Computes a smoothed Dirac delta function based on the signed distance φ.
% Used in level-set methods for interface-localized integrals.
%
% Inputs:
%   phi - signed distance function (vector or matrix)
%   dx  - grid spacing (used to define smoothing width)
%
% Output:
%   res - smoothed delta function values (nonzero only near φ = 0)

    epsilon = 1.5 * dx;        % Width of smoothing band

    res = zeros(size(phi));    % Preallocate output

    % --- Region: φ < -ε (far inside) → delta = 0
    idx = phi < -epsilon;
    res(idx) = 0;

    % --- Region: |φ| ≤ ε (transition zone)
    idx = (-epsilon < phi) & (phi < epsilon);
    res(idx) = (1 / (2 * epsilon)) * (1 + cos(pi * phi(idx) / epsilon));

    % --- Region: φ > ε (far outside) → delta = 0
    idx = phi > epsilon;
    res(idx) = 0;
end
