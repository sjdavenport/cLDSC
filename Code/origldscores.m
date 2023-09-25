function [ ldscores_adjusted, ldscores, chi2 ] = origldscores( n, m, h2, FWHM, mixratio )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'FWHM', 'var' )
   % Default value
   FWHM = 0;
end

if ~exist( 'mixratio', 'var' )
   % Default value
   mixratio = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
if FWHM == 0
    X = randn(n,m);
else
    nsmooth = floor(mixratio*n);
    lat_data = wfield(nsmooth,m);
    smooth_data = convfield(lat_data, FWHM);
    X = [randn(n-nsmooth,m);smooth_data.field];
end

X = X - mean(X);
X = X./std(X,0,1);

% iid beta
beta = ((h2/m)^(1/2))*randn(m,1); % similar to gwash sims, once you normalize
% Email David about the right form of the write it down

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
    loader(j,m, 'Computing ldscores:')
    ldscores(j) = (1/n^2)*X(:,j)'*XXT*X(:,j);
end

ldscores_adjusted = ldscores - (m-ldscores)/(n-2);


end

