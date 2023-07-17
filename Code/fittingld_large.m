%%
h2 = 0.5;
n = 10000;
m = 50000;
% m = 100;

X = randn(n,m);
% How to generate a realistic X for gwas containing 0,1,2. How to get the
% probabilities??
% X = wfield(n,m,'T', 3).field;
% Make binomial random matrix
% X = (randn(n,m) > 0.25);
 
% % X with 0, 1, 2
% X = (randn(n,m) > 0.25) + (randn(n,m) > 0.25);
% % X = smooth_data.field;

X = X - mean(X);
X = X./std(X,0,1);

% lat_data = wfield(n,m);
% FWHM = 3;
% smooth_data = convfield(lat_data, FWHM);
% X = smooth_data.field;
% X = X./std(X,0,1);

% iid beta
beta = ((h2/m)^(1/2))*randn(m,1); % similar to gwash sims, once you normalize

% Actually seems to improve if you makes the betas smooth!
% FWHM = 30;
% smooth_prevarscaled_betas = convfield(wfield(m,1), FWHM);
% global PIloc
% load([PIloc,'Variance/storevars'], 'allvars')
% beta = (smooth_prevarscaled_betas.field/sqrt(allvars(FWHM)))*((h2/m)^(1/2));

% beta = randn(m,1);
e = ((1-h2)^(1/2))*randn(n,1);

phi = X*beta + e;

betahat = zeros(m,1);
for j = 1:m
    betahat(j) = X(:,j)'*phi/n; %/norm(X(:,j))^2 consider normalizing the Xs
end

chi2 = n*betahat.^2;

gwashmn = (mean(chi2) - 1)*(m/n)

%% Other sample 
% Setting this to be 100 whilst setting the original n, m  to be 2000 is an
% interesting example. Seems to indicate that n2 needs to be used.
n2 = 1000;
X = randn(n2,m);
beta = ((h2/m)^(1/2))*randn(m,1); % similar to gwash sims, once you normalize

% beta = randn(m,1);
e = ((1-h2)^(1/2))*randn(n2,1);

phi = X*beta + e;

betahat = zeros(m,1);
for j = 1:m
    betahat(j) = X(:,j)'*phi/n2; %/norm(X(:,j))^2 consider normalizing the Xs
end

chi2_othersample = n2*betahat.^2;

% %%
% betahat(j) = X(:,j)'*phi/n/norm(X(:,j))^2;

XXT = X*X';

ldscores_othersample = zeros(1,m);
for j = 1:m
    loader(j,m, 'Total progress:')
    ldscores_othersample(j) = (1/n2^2)*X(:,j)'*XXT*X(:,j);
end

ldscores_othersample_adjusted = ldscores_othersample - (m-ldscores_othersample)/(n2-2);
ldscores_othersample_adjusted2 = ldscores_othersample_adjusted - (m-ldscores_othersample_adjusted)/(n2-2);

% Unconstrained regression
design = [(ldscores_othersample_adjusted)'*(n2/m), ones(m,1)];
ldsc = (design'*design)^(-1)*design'*chi2_othersample;
ldsc = ldsc(1);

% LD score regression, intercept = 1
design = [(ldscores_othersample_adjusted)'*(n2/m)];
ldsc1 = (design'*design)^(-1)*design'*(chi2_othersample-1);

% Conditional
design = [(ldscores_othersample*(n2/m)-1)'];
(design'*design)^(-1)*design'*(chi2_othersample-1);

% GWASH
gwash = (mean(chi2_othersample) - 1)/mean((ldscores_othersample*(n2/m)-1)')';

gwashmn = (mean(chi2_othersample) - 1)*(m/n2); %This only seems to be correct for Gaussian X, it breaks quite

fprintf('True  | Full LDSC | LDSC intercept 1 |  GWASH   | GWASH m/n\n')
fprintf('%.2f  |   %.2f    |      %.2f        |  %.2f    | %.2f \n', h2, ldsc, ldsc1, gwash, gwashmn)

%% Combining the samples
% Unconstrained LD score regression
design = [(ldscores_othersample_adjusted)'*(n/m), ones(m,1)];

ldsc = (design'*design)^(-1)*design'*chi2;
ldsc = ldsc(1);

% LDSC with the intercept set to 1
design = [(ldscores_othersample_adjusted)'*(n/m)];

ldsc1 = (design'*design)^(-1)*design'*(chi2-1);

% GWASH
gwash = (mean(chi2) - 1)/mean((ldscores_othersample_adjusted*(n/m))');

% GWASH m/n
gwashmn = (mean(chi2) - 1)*(m/n); %This only seems to be correct for Gaussian X, it breaks quite

fprintf('Estimates using the other sample\n')

fprintf('True  | Full LDSC | LDSC intercept 1 |  GWASH   | GWASH m/n\n')
fprintf('%.2f  |   %.2f    |      %.2f        |  %.2f    | %.2f \n', h2, ldsc, ldsc1, gwash, gwashmn)
