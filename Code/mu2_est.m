function [mu1, mu2, S_tilde_mask] = mu2_est(n, m, S_tilde, mu23_version, my_mask, I2)
% mu2 estimation Eq (22) or Eq (27)
mu1 = 1; % mu1 is always 1
if mu23_version == 1
    % Eq(22), full sample
    mu2 = sum(S_tilde.^2) / m - (m - 1) / (n - 1);
else
    % Eq(27), mu.2, I2, subset
    S_tilde_mask = S_tilde .* my_mask;
    mu2 = sum(S_tilde_mask.^2) / m - (mu1^2 * I2 / m) / (n - 1);
end
end