%%
h2 = 0.4;
n = 1000;
m = 2000;

FWHM = 50; mratio = 0.5;
[ ldscores_adjusted, ldscores, chi2 ] = origldscores( n, m, h2, FWHM, mratio );
[ ldsc_full, ldsc_intercept1, ldsc_conditional, gwash, gwashmn] = ...
                                h2ests( n, m, ldscores_adjusted, ldscores, chi2 );
% add the weighted esimator in here!
fprintf('\n')
fprintf('True  | Full LDSC    | LDSC intercept 1  | cLDSC   | GWASH     | GWASH m/n\n')
if ldsc_full(1) > 0
    fprintf('%.2f  |   %.4f     |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc_full(1), ldsc_intercept1, ldsc_conditional, gwash, gwashmn)
else
    fprintf('%.2f  |   %.4f    |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc_full(1), ldsc_intercept1, ldsc_conditional, gwash, gwashmn)  
end
fprintf('LDSC intercept: %.2f \n' ,ldsc_full(2))
% add the weighted esimator in here!

%% Other sample
n2 = 500;
FWHM = 5;
mratio2 = 0.65;
[ ldscores_othersample_adjusted, ldscores_othersample, chi2_othersample ] = ...
                                    origldscores( n2, m, h2, FWHM, mratio2 );
[ ldsc_full, ldsc_intercept1, ldsc_conditional, gwash, gwashmn] = ...
           h2ests( n2, m, ldscores_othersample_adjusted, ldscores_othersample, chi2_othersample );

fprintf('\n')
fprintf('True  | Full LDSC    | LDSC intercept 1  | cLDSC   | GWASH     | GWASH m/n\n')
if ldsc_full(1) > 0
    fprintf('%.2f  |   %.4f     |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc_full(1), ldsc_intercept1, ldsc_conditional, gwash, gwashmn)
else
    fprintf('%.2f  |   %.4f    |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc_full(1), ldsc_intercept1, ldsc_conditional, gwash, gwashmn)  
end
fprintf('LDSC intercept: %.2f \n' ,ldsc_full(2))

%% Using original chi2
[ ldsc_full, ldsc_intercept1, ldsc_conditional, gwash, gwashmn] = ...
           h2ests( n, m, ldscores_othersample_adjusted, ldscores_othersample, chi2 );

fprintf('\n')
fprintf('True  | Full LDSC    | LDSC intercept 1  | cLDSC   | GWASH     | GWASH m/n\n')
if ldsc_full(1) > 0
    fprintf('%.2f  |   %.4f     |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc_full(1), ldsc_intercept1, ldsc_conditional, gwash, gwashmn)
else
    fprintf('%.2f  |   %.4f    |      %.4f       |  %.4f  | %.4f    | %.4f \n', h2, ldsc_full(1), ldsc_intercept1, ldsc_conditional, gwash, gwashmn)  
end
fprintf('LDSC intercept: %.2f \n' ,ldsc_full(2))

%% Other sample 
% Setting this to be 100 whilst setting the original n, m  to be 2000 is an
% interesting example. Seems to indicate that n2 needs to be used.
n2 = 500
X = randn(n2,m);
% lat_data = wfield(n2,m);
% smooth_data = convfield(lat_data, FWHM);
% X = smooth_data.field;

X = X - mean(X);
X = X./std(X,0,1);

% Need to decide whether to recalculate beta here again below!!
% beta = ((h2/m)^(1/2))*randn(m,1); % similar to gwash sims, once you normalize

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

fprintf('True  | Full LDSC | LDSC intercept 1 |  cLDSC  | GWASH   | GWASH m/n\n')
fprintf('%.4f  |   %.4f    |      %.4f        |  %.4f   |  %.4f   | %.4f \n', h2, ldsc, ldsc1, cldsc, gwash, gwashmn)

%% Estimates for the original sample based on these ldscores
plot(ldscores_othersample, ldscores, '*')

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

%% Conditional estimates which are crap
design = [(ldscores_othersample_adjusted*(n/m)-1)'];

(design'*design)^(-1)*design'*(chi2-1)

design = [(ldscores_othersample_adjusted2*(n/m)-1)'];

(design'*design)^(-1)*design'*(chi2-1)

%% Everything below is deprecated
%% Using n
design = [(ldscores_othersample)'*(n/m), ones(m,1)];

(design'*design)^(-1)*design'*chi2

design = [(ldscores_othersample)'*(n/m)];

(design'*design)^(-1)*design'*(chi2-1)

design = [(ldscores_othersample*(n/m)-1)'];

(design'*design)^(-1)*design'*(chi2-1)

%% Estimates using the other sample's ld scores (using n2)
% n2 artificially made things work in certain situations I think - usually its pretty bad!
% there's no motivation for using it i don't think.

design = [(ldscores_othersample)'*(n2/m), ones(m,1)];

(design'*design)^(-1)*design'*chi2

design = [(ldscores_othersample)'*(n2/m)];

(design'*design)^(-1)*design'*(chi2-1)

design = [(ldscores_othersample*(n2/m)-1)'];

(design'*design)^(-1)*design'*(chi2-1)

%%
h2 = 0.4;
nsubj = 10000;
ldscores = randi(250, nsubj, 1)*(1/100);

Y = h2*(ldscores - 1) + 1 + 10*randn(nsubj, 1);

design = [ldscores, ones(nsubj,1)];

(design'*design)^(-1)*design'*Y

%%
% m = 100;

% X = randn(n,m);
% How to generate a realistic X for gwas containing 0,1,2. How to get the
% probabilities??
% X = wfield(n,m,'T', 3).field;
% Make binomial random matrix
% X = (randn(n,m) > 0.45);
 
% X with 0, 1, 2
% X = (randn(n,m) > 0.45) + (randn(n,m) > 0.45);
% X = smooth_data.field;