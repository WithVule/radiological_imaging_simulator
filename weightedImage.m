function slika = weightedImage(ix, iy, mask, E_ph, E, eps_curve, n_d)
    if isscalar(eps_curve)
        w = eps_curve * ones(sum(mask), 1);
    else
        w = interp1(E, eps_curve, E_ph(mask), 'linear', 'extrap');
        w = max(w, 0);
    end

    lin_idx = sub2ind([n_d, n_d], iy(mask), ix(mask));
    slika = accumarray(lin_idx, w, [n_d*n_d, 1]);
    slika = reshape(slika, n_d, n_d);
    slika = flip(slika, 1);
end