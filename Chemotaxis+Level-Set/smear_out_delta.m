function res = smear_out_delta(phi,dx)

    epsilon = 1.5*dx;

    res = zeros(size(phi));

    idx = phi < -epsilon;
    res(idx) = 0;


    idx = (-epsilon < phi) & (phi < epsilon);
    res(idx) = 1/(2*epsilon)+1/(2*phi(idx))*cos(pi*phi(idx)/epsilon);

    idx = phi > epsilon;
    res(idx) = 0;


end 