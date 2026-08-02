function E_ph = photonEnergy(xi, E, N0)

E  = E(:);
N0 = N0(:);
    
f = N0 / sum(N0);
F = cumsum(f);

[F_u, idx] = unique(F, 'stable');
E_u = E(idx);

if F_u(1) > 0
    F_u = [0; F_u];
    E_u = [E_u(1); E_u];
end
    
E_ph = interp1(F_u, E_u, xi, 'linear', 'extrap');

end