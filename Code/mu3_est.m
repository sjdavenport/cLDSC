function mu3 = mu3_est(n, m, S_tilde, mu1, mu2, mu23_version, I2, I3)
% mu3 estimation Eq (30) or Eq (31)
    if mu23_version == 1
        % Eq (30), full sample
        S2 = S_tilde * S_tilde';
        mu3 = sum(sum(S2 * S_tilde)) / m - 3 * (m - 1) * mu2 / (n - 1) - (m - 1) * (m - 2) / (n - 1)^2;
    else
        % Eq (31): mu.3, I3, subset
        mu3 = sum(diag(S_tilde * S_tilde * S_tilde)) / m - 3 * mu2 * (I2 / m) / (n - 1) - (I3 / m) / (n - 1)^2;
    end
end