%%
h2 = 0.2;
n = 500;
m = 1000;

FWHM = 0; mratio = 0.5;
ldsc_full_store = zeros(1,100);
nsim = 1000;
rho = 0.2;
for I = 1:1000
    I
    [ ldscores_adjusted, ldscores, chi2 ] = origldscores( n, m, h2, rho, mratio );
    ldsc_full = h2ests( n, m, ldscores_adjusted, ldscores, chi2 );
    ldsc_full_store(I) =  ldsc_full(1);
end

%%
mean(ldsc_full_store)
std(ldsc_full_store)/sqrt(nsim)