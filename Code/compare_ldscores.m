h2 = 0.4;
n = 1000;
m = 2000;

FWHM = 50; mratio = 1;

[ ldscores_adjusted_1, ldscores_1, chi2_1 ] = origldscores( n, m, h2, FWHM, mratio );
% [ ldscores_adjusted_2, ldscores_2, chi2_2 ] = origldscores( n, m, h2, FWHM, mratio );

%%
histogram(ldscores_1)
mean(ldscores_1)

%%
plot(ldscores_adjusted_1, ldscores_adjusted_2, 'o')

%%
plot(ldscores_1, ldscores_2, 'o')