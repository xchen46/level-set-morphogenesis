function [u,v] = initalize_two_blob(dom_x,dom_y,cx1,cx2,cy1,cy2)


    sigma_u = 0.7;
    
    % Create two Gaussians
    u1 = exp(-((dom_x - cx1).^2 + (dom_y - cy1).^2) / (2*sigma_u^2));
    u2 = exp(-((dom_x - cx2).^2 + (dom_y - cy2).^2) / (2*sigma_u^2));
    u = u1 + u2;
    
    % Initialize v as a smoothed version of u
    sigma_v = 1.5;
    kernel_size = 7;
    [xk, yk] = meshgrid(-kernel_size:kernel_size, -kernel_size:kernel_size);
    G = exp(-(xk.^2 + yk.^2) / (2 * sigma_v^2));
    G = G / sum(G(:));
    v = conv2(u, G, 'same')/10;
    
    u = u(:);
    v = v(:);

end