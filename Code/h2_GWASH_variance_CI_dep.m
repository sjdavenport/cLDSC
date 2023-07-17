function [variance, upper_h2, lower_h2] = h2_GWASH_variance_CI_dep(n, m, h2, mu1, mu2, mu3)
% Asymptotic 95% CI and variance of GWASH h^2 estimate
psi2 = 2 * (m * mu1^2 / (n * mu2) + 2 * mu1 * mu3 * h2 / mu2^2 - h2^2);
upper_h2 = h2 + 1.96 * sqrt(psi2 / n);
lower_h2 = h2 - 1.96 * sqrt(psi2 / n);
variance = psi2;
end