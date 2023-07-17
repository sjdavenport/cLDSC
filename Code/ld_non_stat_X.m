%%
h2 = 0.4;
n = 2000;
m = 2000;
% m = 100;

X1 = randn(n,m/2);
lat_data = wfield(n,m/2);
FWHM = 25;
smooth_data = convfield(lat_data, FWHM);
X2 = smooth_data.field;

X = [X1, X2];
X = X - mean(X);
X = X./std(X,0,1);


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

% %%
% betahat(j) = X(:,j)'*phi/n/norm(X(:,j))^2;

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

cldsc = (design'*design)^(-1)*design'*(chi2-1);

% GWASH (up to O(1/n)) see GWASH supplementary!
gwash = (mean(chi2) - 1)/mean((ldscores_adjusted*(n/m))');

% condtional GWASH
% (mean(chi2) - 1)/mean((ldscores*(n/m)-1)')

gwashmn = (mean(chi2) - 1)*(m/n);

fprintf('True  | Full LDSC    | LDSC intercept 1  | cLDSC   | GWASH     | GWASH m/n\n')
if ldsc > 0
    fprintf('%.2f  |   %.4f     |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc, ldsc1, cldsc, gwash, gwashmn)
else
    fprintf('%.2f  |   %.4f    |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc, ldsc1, cldsc, gwash, gwashmn)  
end