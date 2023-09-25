h2 = 0.5;
n = 2000;
m = 2000;
s
sigma_s_2 = 0.2;

X = randn(n,m);
X = X - mean(X);
X = X./std(X,0,1);
beta = ((h2/m)^(1/2))*randn(m,1); % similar to gwash sims, once you normalize

% beta = randn(m,1);
if (1-h2-sigma_s_2) < 0
    error('1-h2-sigma_s_2 must be greater than 0!')
end
e = ((1-h2-sigma_s_2)^(1/2))*randn(n,1);

phi = X*beta + e;
phi(1:(n/2)) = phi(1:(n/2)) + sqrt(sigma_s_2)/2;
phi((n/2 + 1):n) = phi((n/2 + 1):n) - sqrt(sigma_s_2)/2;

betahat = zeros(m,1);
for j = 1:m
    betahat(j) = X(:,j)'*phi/n; %/norm(X(:,j))^2 consider normalizing the Xs
end

chi2 = n*betahat.^2;

XXT = X*X';

ldscores = zeros(1,m);
for j = 1:m
    loader(j,m, 'Total progress:')
    ldscores(j) = (1/n^2)*X(:,j)'*XXT*X(:,j);
end

ldscores_adjusted = ldscores - (m-ldscores)/(n-2);

% Unconstrained ld score regression
design = [(ldscores_adjusted)'*(n/m), ones(m,1)];

ldsc = (design'*design)^(-1)*design'*chi2;
ldsc = ldsc(1);

% Ld score regression with the intercept set to 1
design = [(ldscores_adjusted)'*(n/m)];

ldsc1 = (design'*design)^(-1)*design'*(chi2-1);

% Conditional
design = [(ldscores*(n/m)-1)'];

(design'*design)^(-1)*design'*(chi2-1);

% GWASH (up to O(1/n)) see GWASH supplementary!
gwash = (mean(chi2) - 1)/mean((ldscores_adjusted*(n/m))');

% condtional GWASH
% (mean(chi2) - 1)/mean((ldscores*(n/m)-1)')

gwashmn = (mean(chi2) - 1)*(m/n);

fprintf('True  | Full LDSC | LDSC intercept 1 |  GWASH   | GWASH m/n\n')
fprintf('%.2f  |   %.2f    |      %.2f        |  %.2f    | %.2f \n', h2, ldsc, ldsc1, gwash, gwashmn)

%%
n2 = 500;
X = randn(n2,m);
X = X - mean(X);
X = X./std(X,0,1);
XXT = X*X';

e = ((1-h2-sigma_s_2)^(1/2))*randn(n2,1);

phi = X*beta + e;
phi(1:(n2/2)) = phi(1:(n2/2)) + sqrt(sigma_s_2)/2;
phi((n2/2 + 1):n2) = phi((n2/2 + 1):n2) - sqrt(sigma_s_2)/2;

betahat = zeros(m,1);
for j = 1:m
    betahat(j) = X(:,j)'*phi/n2; %/norm(X(:,j))^2 consider normalizing the Xs
end

chi2_othersample = n2*betahat.^2;

ldscores_othersample = zeros(1,m);
for j = 1:m
    loader(j,m, 'Total progress:')
    ldscores_othersample(j) = (1/n2^2)*X(:,j)'*XXT*X(:,j);
end

ldscores_othersample_adjusted = ldscores_othersample - (m-ldscores_othersample)/(n2-2);

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

%% Using n
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
fprintf('%.4f  |   %.4f    |      %.4f        |  %.4f    | %.4f \n', h2, ldsc, ldsc1, gwash, gwashmn)


%% Large dataset
h2 = 0.2;
n = 2000;
m = 2000;

sigma_s_2 = 0.2;

X = randn(n,m);
X = X - mean(X);
X = X./std(X,0,1);
beta = ((h2/m)^(1/2))*randn(m,1); % similar to gwash sims, once you normalize

% beta = randn(m,1);
if (1-h2-sigma_s_2) < 0
    error('1-h2-sigma_s_2 must be greater than 0!')
end
e = ((1-h2-sigma_s_2)^(1/2))*randn(n,1);

phi = X*beta + e;
phi(1:(n/2)) = phi(1:(n/2)) + sqrt(sigma_s_2)/2;
phi((n/2 + 1):n) = phi((n/2 + 1):n) - sqrt(sigma_s_2)/2;

betahat = zeros(m,1);
for j = 1:m
    betahat(j) = X(:,j)'*phi/n; %/norm(X(:,j))^2 consider normalizing the Xs
end

chi2 = n*betahat.^2;

gwashmn = (mean(chi2) - 1)*(m/n)