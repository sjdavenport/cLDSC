% Information needed for simulating dataset
h2 = 0.2; % Heritability
m = 250; % number of SNP
n = 100; % number of subject
nsim = 50; % number of dataset simulated, if nsim=1, only 1 dataset is generated

b.dist = 2; % 1: beta is from N(0,1); 2: beta is the mixture of N(0,1) and 0's
X.norm = false; % X.norm=true if X is normal distribution; X.norm=false if X is binomial distribution.
rho = 0.2; % rho>0 is the AR correlation coefficient

if (~X.norm) % if X is binomial,
    p = 0.1;
    q = 1 - p; % binomial success and failure rates, could be other values;
    rep.num = 10; % for obtaining covariance matrix;
end

if (b.dist == 2) % if beta is mixed
    null.b.prop = 0.9; % 90% of beta are 0, the other 10% beta are from N(0,1)
end

FWHM = 3;
for I = 1:rep.num
    lat_data = wfield(true(n,1), m);
    X = convfield(lat_data, FWHM);
    correlated_binom = X > p;
end
SIGMA = covmate.sum / rep.num; % SIGMA is an average of rep.num

% Generate variance-covariance matrix for X, assuming AR
if (X.norm) % X is normal distribution
    Sigma = rho .^ abs(outer(1:m, 1:m, '-')); % correlation matrix
    SIGMA = sqrt(diag(1:m)) * Sigma * sqrt(diag(1:m)); % fake a covariance matrix
else % X is binomial
    covmate.sum = zeros(m);
    for i = 1:rep.num
        normal.cop = copularnd('Gaussian', rho, m);
        u = normcdf(normal.cop);
        correlated.binom = icdf('Binomial', u, 2, p);
        covmate.sum = covmate.sum + covmate(correlated.binom);
    end
    SIGMA = covmate.sum / rep.num; % SIGMA is an average of rep.num
end

% Generate beta
if (b.dist == 1) % beta is from N(0,1)
    b_star = randn(m, 1);
else % beta is mixture
    b_star = zeros(m, 1);
    null.b.indx = sort(randsample(1:m, round(m * null.b.prop), false));
    b_star(setdiff(1:m, null.b.indx)) = randn(round((1 - null.b.prop) * m), 1);
end

% Normalize beta
b = b_star * sqrt(h2) / sqrt(b_star' * SIGMA * b_star);
sigma2.eps = 1 - h2;

% Prepare for subset for computational efficiency,
% when using the full dataset, this isn't needed;
% n.subdiag is equivalent to q in Table 1
mask = abs(outer(1:m, 1:m, '-')) <= n.subdiag;
I2 = sum(mask(:)) - m;
array.mask = mask;

[ind1, ind2, ind3] = ndgrid(1:m, 1:m, 1:m);
ind = [ind1(:), ind2(:), ind3(:)];
valid.ind = (ind(:, 1) ~= ind(:, 2)) & (ind(:, 1) ~= ind(:, 3)) & (ind(:, 2) ~= ind(:, 3))

%%
% Create variables to hold the results
z_mat = zeros(m, nsim);

% For all the following:
% v1: full sample; v2: subset

mu2_v1 = zeros(1, nsim);
mu2_v2 = zeros(1, nsim);
mu3_v1 = zeros(1, nsim);
mu3_v2 = zeros(1, nsim);

% GAWSH h^2
h2_GWASH_v1 = zeros(1, nsim);
h2_GWASH_v2 = zeros(1, nsim);

% Asymptotic 95% CI and variance for GAWSH h^2
h2_var_v1 = zeros(1, nsim);
h2_upper_v1 = zeros(1, nsim);
h2_lower_v1 = zeros(1, nsim);

h2_var_v2 = zeros(1, nsim);
h2_upper_v2 = zeros(1, nsim);
h2_lower_v2 = zeros(1, nsim);

% Generate totally nsim dataset and obtain parameter estimates from each dataset
for i = 1:nsim
    % Generate data
    if X_norm
        X = mvnrnd(zeros(1, m), SIGMA, n);
    else
        normal_cop = copularnd('Normal', rho, m);
        u = copulainv('Binormal', normal_cop, p);
        X = icdf('Binomial', u, 2, p);
    end
    
    eps = normrnd(0, sqrt(sigma2_eps), n, 1);
    y = X * b + eps;
    
    X_s = zscore(X);
    y_s = zscore(y);
    
    z_mat(:, i) = X_s' * y_s / sqrt(n - 1); % Eq (16) in vector form
    
    % Parameter estimates
    S = (X_s' * X_s) / (n - 1); % Eq (17), sample covariance matrix
    
    % version 1: Full data
    m_obj = mu2_est(n, m, S, 1, 0, 0);
    mu1 = m_obj.mu1; % always 1
    mu2_v1(i) = m_obj.mu2;
    
    % GAWSH estimate h^2 under full data
    h2_GWASH_v1(i) = h2_GWASH_est_dep(z_mat(:, i), mu1, mu2_v1(i), n, m);
    
    mu3_v1(i) = mu3_est(n, m, S, mu1, mu2_v1(i), 1, 0, 0).mu3;
    
    % SE estimate for GAWSH estimate h^2 under entire data
    est_obj = h2_GWASH_variance_CI_dep(n, m, h2_GWASH_v1(i), mu1, mu2_v1(i), mu3_v1(i));
    h2_var_v1(i) = est_obj.variance;
    h2_upper_v1(i) = est_obj.upper_h2;
    h2_lower_v1(i) = est_obj.lower_h2;
    
    % version 2: subset
    m_obj = mu2_est(n, m, S, 2, array_mask, I2);
    mu2_v2(i) = m_obj.mu2;
    S_mask = m_obj.S_tilde_mask;
    
    % GAWSH estimate h^2 under subset
    h2_GWASH_v2(i) = h2_GWASH_est_dep(z_mat(:, i), mu1, mu2_v2(i), n, m);
    
    mu3_v2(i) = mu3_est(n, m, S_mask, mu1, mu2_v2(i), 2, I2, I3).mu3;
    
    % SE estimate for GAWSH estimate h^2 under subset
    est_obj = h2_GWASH_variance_CI_dep(n, m, h2_GWASH_v2(i), mu1, mu2_v2(i), mu3_v2(i));
    h2_var_v2(i) = est_obj.variance;
    h2_upper_v2(i) = est_obj.upper_h2;
    h2_lower_v2(i) = est_obj.lower_h2;
    
end

%%
% Report results

% Based on full sample, v1
disp(['h2.GWASH.v1: ', num2str(mean(h2_GWASH_v1))]);
disp(['Empirical S.E of GWASH h^2: ', num2str(std(h2_GWASH_v1))]);

% mu2
disp(['mu2.v1: ', num2str(mean(mu2_v1))]);
disp(['Empirical S.E of mu2: ', num2str(std(mu2_v1))]);

% mu3
disp(['mu3.v1: ', num2str(mean(mu3_v1))]);
disp(['Empirical S.E of mu3: ', num2str(std(mu3_v1))]);

% Asymptotic variance and 95% CI
disp(['asymptotic variance of GWASH h^2: ', num2str(mean(h2_var_v1))]);
disp(['asymptotic S.E of GWASH h^2: ', num2str(mean(sqrt(h2_var_v1) / n))]);
disp(['lower: ', num2str(mean(h2_lower_v1(~isnan(h2_lower_v1))))]);
disp(['upper: ', num2str(mean(h2_upper_v1(~isnan(h2_upper_v1))))]);

% Based on subset, v2
disp(['h2.GWASH.v2: ', num2str(mean(h2_GWASH_v2))]);
disp(['Empirical S.E of GWASH h^2: ', num2str(std(h2_GWASH_v2))]);

% mu2
disp(['mu2.v2: ', num2str(mean(mu2_v2))]);
disp(['Empirical S.E of mu2: ', num2str(std(mu2_v2))]);

% mu3
disp(['mu3.v2: ', num2str(mean(mu3_v2))]);
disp(['Empirical S.E of mu3: ', num2str(std(mu3_v2))]);

% Asymptotic variance and 95% CI
disp(['asymptotic variance of GWASH h^2: ', num2str(mean(h2_var_v2))]);
disp(['asymptotic S.E of GWASH h^2: ', num2str(mean(sqrt(h2_var_v2) / n))]);
disp(['lower: ', num2str(mean(h2_lower_v2(~isnan(h2_lower_v2))))]);
disp(['upper: ', num2str(mean(h2_upper_v2(~isnan(h2_upper_v2))))]);
