function q_out = quatmult(p, q)

eta_p = p(4);
eps_p = p(1:3);
eta_q = q(4);
eps_q = q(1:3);

q_out = [eta_p*eps_q + eta_q*eps_p + (vect2cross(eps_p)*eps_q); (eta_p*eta_q - (eps_p'*eps_q))];

end

