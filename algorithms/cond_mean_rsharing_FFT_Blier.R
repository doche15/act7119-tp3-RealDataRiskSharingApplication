### ACT-7119: Modèles de risque avec dépendance et mesures de risque
### TP3 : Partage de risque
### Calcul de la règle de l'espérance conditionnelle

### Code tiré de l'annexe A2 de Blier-Wong et al. (2025)

n_participants <- 24500
kmax <- 2^14
lam <- list()
fc <- list()
mu <- list()
lambdas <- params_data$lambda[1:n_participants]

h <- 250
# Assign parameters
for(i in 1:n_participants) {
    lam[[i]] <- lambdas[i]
    fc[[i]] <- discretize_gamma_dist(kmax, h, params_data$alpha[i],
                                     params_data$beta[i])
    print(i)
}
dft_fx <- list()
phic <- list()
cm <- list()
for(i in 1:n_participants) {
    dft_fx[[i]] <- exp(lam[[i]] * (fft(fc[[i]])- 1))
    phic[[i]] <- fft(c(1:(kmax-1) * h * fc[[i]][-1], 0))
    print(i)
}
dft_fs <- Reduce("*", dft_fx)
fs <- Re(fft(dft_fs, inverse = TRUE))/kmax
e1 <- exp(-2i*pi*(0:(kmax-1))/kmax)
for(i in 1:n_participants) {
    dft_mu <- e1 * phic[[i]] * lam[[i]] * dft_fs
    mu[[i]] <- Re(fft(dft_mu, inverse = TRUE))/kmax
    cm[[i]] <- mu[[i]]/fs
    print(i)
}
cm_tot <- Reduce("+", cm)
## Pour que la validation marche, il faut visualiser où il y a beaucoup de
## masse dans fs avec
plot(0:(kmax-1) * h, fs)

cbind(cm_tot[1 + seq(2.8e6, 3.2e6, 100) / h], # Validation
(seq(2.8e6, 3.2e6, 100)))


purepremium <- params_data$lambda * params_data$alpha / params_data$beta

xx <- 9
cm[xx]
purepremium[xx] # Écart vu le pas de discrétisation



fss <- fs[(cm[[9]] > 0) & (cm[[9]] <1000)]
cmm <- cm[[9]][(cm[[9]] > 0) & (cm[[9]] <1000)]
plot(sort(cmm), cumsum(fss[order(cmm)]))
abline(v = purepremium[9], col = "red")
