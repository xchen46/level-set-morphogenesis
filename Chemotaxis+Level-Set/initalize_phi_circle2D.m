%% Initialized phi as a 2D circle

function phi = initalize_phi_circle2D(dom_x,dom_y,cx,cy,r)

    phi = sqrt((dom_x-cx).^2+(dom_y-cy).^2)-r;
    phi = phi(:);

end
