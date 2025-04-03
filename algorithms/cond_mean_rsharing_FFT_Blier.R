### ACT-7119: Modèles de risque avec dépendance et mesures de risque
### TP3 : Partage de risque
### Calcul de la règle de l'espérance conditionnelle

### Code tiré et adapté de l'annexe A2 de Blier-Wong et al. (2025)

library(actuar)
## Discrétisation

params_data <- read.csv("params_belgian.csv")

discrete_sev <- function(nfft, h, method, alpha, beta)
{
    f_sev <- discretize(pgamma(x, alpha, beta),
                        method = method,
                        from = 0, to = (nfft - 1) * h, step = h,
                        lev = levgamma(x, alpha, beta))

    if (method == "upper")
        return(c(f_sev, 0)) # Pour avoir une longueur de nfft (bug d'actuar)

    f_sev
}

## Donne le même résultat que notre fonction maison
# cbind(discrete_sev(kmax, h, "upper", params_data$alpha[2], params_data$beta[2]),
#       discretize_gamma_dist(kmax, h, params_data$alpha[2], params_data$beta[2],
#                             m = "upper"))
# cbind(discrete_sev(kmax, h, "lower", params_data$alpha[2], params_data$beta[2]),
#       discretize_gamma_dist(kmax, h, params_data$alpha[2], params_data$beta[2],
#                             m = "lower"))

## Mean-preserving
h <- 25
kmax <- 2^15
ff <- discrete_sev(kmax, h, "unbiased", params_data$alpha[2], params_data$beta[2])
sum(0:(kmax - 1) * h * ff)
params_data$alpha[2] / params_data$beta[2]

compute_cm <- function(kmax, h, params_data, discr_method)
{
    # kmax: doit être 2^m
    # h: pas de discrétisation
    # params_data: df de paramètres avec colonnes $alpha, $beta et $lambda.

    n_participants <- nrow(params_data)
    lam <- list()
    fc <- list()
    mu <- list()

    total_start <- Sys.time()
    # Assign parameters
    debut <- Sys.time()
    for(i in 1:n_participants) {
        lam[[i]] <- lambdas[i]
        fc[[i]] <- discrete_sev(kmax, h, discr_method, params_data$alpha[i],
                                params_data$beta[i])
        if (i %% 100 == 0)
            print(paste0("Boucle 1 (discretisation) :",
                         round(i / n_participants * 100, 2),
                         " % accompli"))
    }
    end <- Sys.time()
    print(paste0("Boucle 1 (discretisation):", round(difftime(end, debut,
                                                              units = "secs"), 1),
                 " secondes."))

    dft_fx <- list()
    phic <- list()
    cm <- list()

    debut <- Sys.time()
    for(i in 1:n_participants) {
        dft_fx[[i]] <- exp(lam[[i]] * (fft(fc[[i]])- 1))
        phic[[i]] <- fft(c(1:(kmax-1) * h * fc[[i]][-1], 0))
        if (i %% 100 == 0)
            print(paste0("Boucle 2 (phi c) :",
                         round(i / n_participants * 100, 2),
                         " % accompli"))
    }

    end <- Sys.time()
    print(paste0("Boucle 2 (phi c):", round(difftime(end, debut,
                                                     units = "secs"), 1),
                 " secondes."))


    dft_fs <- Reduce("*", dft_fx)
    fs <- Re(fft(dft_fs, inverse = TRUE))/kmax
    e1 <- exp(-2i*pi*(0:(kmax-1))/kmax)

    debut <- Sys.time()
    for(i in 1:n_participants) {
        dft_mu <- e1 * phic[[i]] * lam[[i]] * dft_fs
        mu[[i]] <- Re(fft(dft_mu, inverse = TRUE))/kmax
        cm[[i]] <- mu[[i]]/fs
        if (i %% 100 == 0)
            print(paste0("Boucle 3 (inversion) :",
                         round(i / n_participants * 100, 2),
                         " % accompli"))
    }



    cm_tot <- Reduce("+", cm)
    end <- Sys.time()

    print(paste0("Boucle 3 (inversion):", round(difftime(end, debut,
                                                         units = "secs"), 1),
                 " secondes."))

    total_end <- Sys.time()

    print(paste0("Temps total: ", round(difftime(total_end, total_start,
                                                 units = "secs"), 1),
                 " secondes."))

    list("contrib" = cm, "tot_contrib" = cm_tot, "fs" = fs)
}

# h <- 25
# kmax <- 2^15
# res <- compute_cm(kmax, h, params_data[1:5000, ], "unbiased")
#
# sum(res$fs * (0:(kmax-1) * h))
#
# sum(params_data$lambda[1:5000] * params_data$alpha[1:5000] /
#         params_data$beta[1:5000])
#
# ## Pour que la validation marche, il faut visualiser où il y a beaucoup de
# ## masse dans fs avec
# plot(0:(kmax-1) * h, res$fs)
#
# cbind(res$tot_contrib[1 + seq(4e5, 6e5, 10) / h], # Validation
# ( seq(4e5, 6e5, 10)))
#
#
# purepremium <- params_data$lambda * params_data$alpha / params_data$beta
#
# xx <- 9
# contribxx <- unlist(res$contrib[xx])
# contribxx[5e5 / h] # Faut regarder la validation de la convergence là où on a bcp de
#             # masse pour fS.
# purepremium[xx]
