function smooth_u = bolb_rand(x,membrane_size)

[X, Y] = meshgrid(x, x);
smooth_u = zeros(size(X));

num_blobs = 6;
% rng(1); % for reproducibility

for n = 1:num_blobs
    x0 = membrane_size/2 * rand()+membrane_size/4;
    y0 = membrane_size/2 * rand()+membrane_size/4;
    amp = 0.05 * rand();
    r = 0.1 + 0.1 * rand();  % radius of blob

    smooth_u = smooth_u + amp * exp(-((X - x0).^2 + (Y - y0).^2) / r^2);
end

smooth_u = 0.1 + smooth_u;
% smooth_u = smooth_u / sum(smooth_u(:));  % normalize

% % --- Apply smoothing filter ---
% smooth_u = conv2(smooth_u, 'same');

end