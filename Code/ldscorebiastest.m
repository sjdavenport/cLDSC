%% Computing the bias of the ldscore estimators
h2 = 0.2;
n = 500;
m = 1000;

nsim = 1000;
ldsc_full_store = zeros(1,nsim);
ldsc_inter1 = zeros(1,nsim);
rho = 0.2;
for I = 1:1000
    I
    [ ldscores_adjusted, ldscores, chi2 ] = origldscores( n, m, h2, rho );
    [ldsc_full, ld_inter1] = h2ests( n, m, ldscores_adjusted, ldscores, chi2 );
    ldsc_full_store(I) =  ldsc_full(1);
    ldsc_inter1(I) =  ld_inter1;
end

%%
mean(ldsc_full_store)
std(ldsc_full_store)/sqrt(nsim)
mean(ldsc_inter1)
std(ldsc_inter1)/sqrt(nsim)