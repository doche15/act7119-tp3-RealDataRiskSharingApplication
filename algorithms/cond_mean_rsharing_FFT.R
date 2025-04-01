###
### Travail 3, ACT-7119
### Algorithme 1 de [Blier-Wong et al., 2025]
###
### [Blier-Wong et al., 2025] Blier-Wong, C., Cossette, H., and Marceau, E.
### (2025). Efficient evaluation of risk allocations. Insurance : Mathematics
### and Economics.
###
##


discretize_gamma_dist <- function(n, h, param1, param2, m = "lower", nn = 27000)
{
    # Référence : Notes de cours du chapitre 1, partie 1, du cours ACT-3000.
    fmp <- numeric(n)
    nn <- min(nn / h, n)
    seq_n <- seq(0, nn - 1)
    vecteur_Frep <- pgamma(seq_n[2:nn]* h, param1, param2)

    fmp[1] <- pgamma(h, param1, param2)
    fmp[2:(nn-2)] <- head(diff(vecteur_Frep), -1)
    fmp[1:(nn - 1)] <- c(fmp[1:(nn-2)], 1 - sum(fmp[1:(nn-2)]))


    if (m == "lower")
    {
        return(c(0, head(fmp, -1)))
    }
    else if (m == "upper")
    {
        return(fmp)
    }
    else
    {
        stop("Méthode non reconnue")
    }
}

ff <- discretize_gamma_dist(2^18, 0.001, 0.1, 0.05)
sum(ff * 0.001 * (0:(2^18 - 1)))
evaluate_pgf_Poisson <- function(t, lambda)
{
    exp( lambda * (t - 1) )
}


#
# nfft        : valeur numérique.
# param_lamda : vecteur.
# param_alpha : vecteur.
# param_beta  : vecteur.
# h           : valeur numérique.
# m           : "lower" ou "upper".
#
calculate_cond_mean_rsharing_FFT <- function(nfft, param_lambda, param_alpha,
                                             param_beta, h, m = "lower")
{
    nfft_vec <- seq(0, nfft - 1)

    fft_fS <- rep(1, nfft)

    phiBi <- matrix(numeric(nfft * length(param_lambda)),
                    nrow = length(param_lambda))
    #browser()
    # Étapes 1-3, complétées différemment. Le but est de retrouver fS.
    for (i in seq_along(param_lambda))
    {
        # Selon les notes de cours du chapitre 1, partie 2, du cours ACT-3000.
        fB <- discretize_gamma_dist(nfft, h, param_alpha[i], param_beta[i],
                                    m = m)
        phiBi[i,] <- fft(h * (nfft_vec + 1) * fB[nfft_vec + 1])
        # print(sum(nfft_vec * h * fB)) # OK
        # print(sum(fB)) # OK
        fft_fB <- fft(fB)
        fft_fX <- evaluate_pgf_Poisson(fft_fB, param_lambda[i])

        fft_fS <- fft_fS * fft_fX
    }

    fS <- Re(fft(fft_fS, inverse = TRUE)) / nfft

    print("Somme de fS :")
    print(sum(fS)) # OK
    print("Espérance de S :")
    print(sum(nfft_vec * h * fS)) # OK
    print("Variance de S :")
    print(sum(((nfft_vec * h)^2) * fS) - (sum(nfft_vec * h * fS)^2))

    # e_vec <- exp(1i * 2 * pi * (nfft_vec / nfft)) # "1i" est i si j'ai bien
    # # compris.
    #
    # cond_mean_rsharing <- matrix(numeric(nfft * length(param_lambda)),
    #                              nrow = length(param_lambda))
    # for (i in seq_along(param_lambda))
    # {
    #     phiB <- (nfft_vec + 1) * dgamma((nfft_vec + 1), param_alpha[i],
    #                                     param_beta[i])
    #
    #     fft_mu <- param_lambda[i] * e_vec * phiB * fft_fS
    #
    #     mu <- Re(fft(fft_mu, inverse = TRUE)) / nfft
    #
    #     cond_mean_rsharing[i, ] <- mu / fS
    # }

    # Essai, déboguage
    e_vec <- exp(-2i * pi * ((nfft_vec) / (nfft)))
    # e_vec <- exp(1i * 2 * pi * (h * nfft_vec))
    # Définition cohérente avec l'aide de la fonction "fft".

    cond_mean_rsharing <- matrix(numeric(nfft * length(param_lambda)),
                                 nrow = length(param_lambda))
    j <- 0
    for (i in seq_along(param_lambda))
    {
        j <- j + 1
        # phiB <- fft(h * (nfft_vec + 1) * dgamma(h * (nfft_vec + 1), param_alpha[i],
        #                                         param_beta[i]))
        # # phiB <- fft(h * nfft_vec * dgamma(h * nfft_vec, param_alpha[i], param_beta[i]))

        fft_mu <- param_lambda[i] * e_vec * phiBi[i,] * fft_fS

        mu <- Re(fft(fft_mu, inverse = TRUE)) / nfft

        cond_mean_rsharing[i, ] <- mu / fS
        print(j)
    }

    cond_mean_rsharing
}


# Tests
params_data <- read.csv("params_belgian.csv")


# Espérance de S
sum(params_data[1:3, 1] * (params_data[1:3, 2] / params_data[1:3, 3]))
sum(params_data$lambda[1:3] *
        ( params_data$alpha[1:3] / (params_data$beta[1:3]^2) ) *
        (1 + params_data$alpha[1:3]))

params_data[1:2, 1] * (params_data[1:2, 2] / params_data[1:2, 3])

nfft <- 2^14
h <- 0.1

res <- calculate_cond_mean_rsharing_FFT(nfft, params_data[1:2, 1],
                                        params_data[1:2, 2],
                                        params_data[1:2, 3], h)

dim(res)
colSums(res[,1:10]) / h
params_data[1:2,]
params_data[1:2, 1] * (params_data[1:2, 2] / params_data[1:2, 3])
t(res)


all.equal(colSums(res)[1:50], (h * seq(0, nfft - 1))[1:50], tolerance = 10^(-4))
# Voir page 126 de [Blier-Wong et al., 2025] pour commentaire sur la portée de
# ce test (test peu bon pour les queues de distribution, car pmf de S très
# faible).

sum(params_data[1:50, 1] * (params_data[1:50, 2] / params_data[1:50, 3]))

sum(params_data$lambda[1:50] *
        ( params_data$alpha[1:50] / (params_data$beta[1:50]^2) ) *
        (1 + params_data$alpha[1:50]))

nfft <- 2^18
h <- 0.1

time1 <- Sys.time()
res <- calculate_cond_mean_rsharing_FFT(nfft, params_data[1:30, 1],
                                        params_data[1:30, 2],
                                        params_data[1:30, 3], h)
time2 <- Sys.time()

time2 - time1

sum(res[, 2])

colSums(res[, 1:10])

# sum(params_data[, 1] * (params_data[, 2] / params_data[, 3]))
