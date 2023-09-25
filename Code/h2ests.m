function [ ldsc_full, ldsc_intercept1, ldsc_conditional, gwash, gwashmn] = ...
                          h2ests( n, m, ldscores_adjusted, ldscores, chi2 )
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
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Unconstrained ld score regression
design = [(ldscores_adjusted)'*(n/m), ones(m,1)];
% 
ldsc_full = (design'*design)^(-1)*design'*chi2;
% ldsc = ldsc(1);

% Ld score regression with the intercept set to 1
design = [(ldscores_adjusted)'*(n/m)];

ldsc_intercept1 = (design'*design)^(-1)*design'*(chi2-1);

% Conditional
design = [(ldscores*(n/m)-1)'];
ldsc_conditional = (design'*design)^(-1)*design'*(chi2-1);

% GWASH (up to O(1/n)) see GWASH supplementary!
gwash = (mean(chi2) - 1)/mean((ldscores_adjusted*(n/m))');

% condtional GWASH
% (mean(chi2) - 1)/mean((ldscores*(n/m)-1)')

gwashmn = (mean(chi2) - 1)*(m/n);

end

