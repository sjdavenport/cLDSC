function h2_GWASH = h2_GWASH_est_dep(u, mu1, mu2, n, m)
% GWASH h^2 estimate, Eq (20)
    s2 = mean(u.^2); % Eq. (21)
    h2_GWASH = (m * mu1 / (n * mu2)) * (s2 - mu1); % mu = 1 always, Eq. 20, GWASH
end
