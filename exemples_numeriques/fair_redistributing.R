###
### Travail 3 du cours ACT-7119
###
### Exemple numérique de la propriété de fair redistributing
###
##

library(tidyr)
library(ggplot2)

source("algorithms/cond_mean_rsharing_FFT_Blier.R")
source("algorithms/mean_prop_rsharing.R")
source("algorithms/uniform_rsharing.R")

summary(params_data)

# trouve les trois risques représentatifs ----
# Importation des données de sous-pools
subpool1 <- read.csv("donnees-sous-pools/sous-pool-1-k3--2025-04-04.csv")
subpool2 <- read.csv("donnees-sous-pools/sous-pool-2-k3--2025-04-04.csv")
subpool3 <- read.csv("donnees-sous-pools/sous-pool-3-k3--2025-04-04.csv")


# Constantes
SUBPOOL1_SIZE <- nrow(subpool1)
SUBPOOL2_SIZE <- nrow(subpool2)
SUBPOOL3_SIZE <- nrow(subpool3)

RISK_ECH_SIZE <- 1000
PURE_PREMIUM_RATE_RISK_CHOICE <- 0.05


# Échantillonnage des 1 000 risques des sous-pools
set.seed(11)

subpool1_ech <- subpool1[sample(seq(SUBPOOL1_SIZE), RISK_ECH_SIZE),]
subpool2_ech <- subpool2[sample(seq(SUBPOOL2_SIZE), RISK_ECH_SIZE),]
subpool3_ech <- subpool3[sample(seq(SUBPOOL3_SIZE), RISK_ECH_SIZE),]

prime_pure_moyennes_ech <-c(mean(subpool1_ech$prime_pure),
                            mean(subpool2_ech$prime_pure),
                            mean(subpool3_ech$prime_pure))


str(subpool1_ech)

# Choix du risque typique analysé dans chaque échantillon de sous-pool
id_risk_ech_subpool1 <- subpool1_ech$id[sample(
  which(
    (1 - PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[1] <= subpool1_ech$prime_pure &
      subpool1_ech$prime_pure <= (1 + PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[1]
  ),
  1
)]

id_risk_ech_subpool2 <- subpool2_ech$id[sample(
  which(
    (1 - PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[2] <= subpool2_ech$prime_pure &
      subpool2_ech$prime_pure <= (1 + PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[2]
  ),
  1
)]

id_risk_ech_subpool3 <- subpool3_ech$id[sample(
  which(
    (1 - PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[3] <= subpool3_ech$prime_pure &
      subpool3_ech$prime_pure <= (1 + PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[3]
  ),
  1
)]

param_indiv <- params_data[c(id_risk_ech_subpool1,
                             id_risk_ech_subpool2,
                             id_risk_ech_subpool3),]

# crée pool initial ----
pool_init <- param_indiv

# convolution des risques 2 et 3 ----
# on forme le pool redistribué tel que (Y1, Y2) = (X1, X2 + X3)
lambda_1 <- param_indiv[1,]$lambda
lambda_2 <- param_indiv[2,]$lambda
lambda_3 <- param_indiv[3,]$lambda

lambda_23 <- lambda_2 + lambda_3

nfft <- kmax
h <- 1

fb2 <- discrete_sev(nfft, h, method = 'unbiased',
                    alpha = param_indiv[2,]$alpha,
                    beta = param_indiv[2,]$beta)
sum(fb2)

fb3 <- discrete_sev(nfft, h, method = 'unbiased',
                    alpha = param_indiv[3,]$alpha,
                    beta = param_indiv[3,]$beta)
sum(fb3)

fc <- (lambda_2 / (lambda_2 + lambda_3)) * fb2 +
  (lambda_3 / (lambda_2 + lambda_3)) * fb3
sum(fc)

phic <- fft(fc)

phix23 <- exp(lambda_23 * (phic - 1))

fx23 <- Re(fft(phix23, inverse = T))/nfft
sum(fx23)

s <- 100

# espérance conditionnelle ----

# NOTE : La fonction suivante est gravement hardcoded par souci de temps, mais ça
# se généralise facilement.
compute_cm_fair_redistributing <- function(kmax, h, params_data, discr_method, fc_combined, lambda_combined)
{
  # kmax: doit être 2^m
  # h: pas de discrétisation
  # params_data: df de paramètres avec colonnes $alpha, $beta et $lambda.
  # fc_combined : fonction de densité de la sévérité des risques combinés
  # lambda_combined : lambda des risques combinés
  
  n_participants <- nrow(params_data)
  lam <- list()
  fc <- list()
  mu <- list()
  
  total_start <- Sys.time()
  # Assign parameters
  debut <- Sys.time()
  for(i in 1:n_participants) {
    lam[[i]] <- params_data$lambda[i]
    fc[[i]] <- discrete_sev(kmax, h, discr_method, params_data$alpha[i],
                            params_data$beta[i])
    if (i %% 100 == 0)
      print(paste0("Boucle 1 (discretisation) :",
                   round(i / n_participants * 100, 2),
                   " % accompli"))
  }
  
  lam[[2]] <- lambda_combined
  fc[[2]] <- fc_combined
  
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
  
  dft_fx[[2]] <- exp(lam[[2]] * (fft(fc[[2]])- 1))
  phic[[2]] <- fft(c(1:(kmax-1) * h * fc[[2]][-1], 0))
  
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
  
  dft_mu <- e1 * phic[[2]] * lam[[2]] * dft_fs
  mu[[2]] <- Re(fft(dft_mu, inverse = TRUE))/kmax
  cm[[2]] <- mu[[2]]/fs
  
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

## calcul des contributions du pool initial ----
cm_initial <- compute_cm(kmax = kmax,
                         h = h,
                         params = pool_init,
                         discr_method = 'unbiased')

## calcul des contributions du pool redistribué ----
cm_redist <- compute_cm_fair_redistributing(kmax = kmax,
                                            h = h,
                                            params = pool_init[1,],
                                            discr_method = 'unbiased',
                                            fc_combined = fc,
                                            lambda_combined = lambda_23)

## verif ----
val_range = seq(2e4)
sum(cm_initial$fs[val_range])
sum(cm_redist$fs[val_range])

fs_initial <- cm_initial$fs[val_range]
sum(fs_initial)
fs_redist <- cm_redist$fs[val_range]
sum(fs_redist)

# contributions du pool initial
contrib_id1_init <- unlist(cm_initial$contrib[1])[val_range]
contrib_id2_init <- unlist(cm_initial$contrib[2])[val_range]
contrib_id3_init <- unlist(cm_initial$contrib[3])[val_range]
# prime pure
sum(contrib_id1_init * fs_initial)
sum(contrib_id2_init * fs_initial)
sum(contrib_id3_init * fs_initial)

contrib_id1_init[s + 1]
contrib_id2_init[s + 1]
contrib_id3_init[s + 1]

# contributions du pool redistribué
contrib_id1_redist <- unlist(cm_redist$contrib[1])[val_range]
contrib_id2_redist <- unlist(cm_redist$contrib[2])[val_range]
# prime pure
sum(contrib_id1_redist * fs_redist)
sum(contrib_id2_redist * fs_redist)
# pour s = 1
contrib_id1_redist[s + 1]
contrib_id2_redist[s + 1]

# prorata ----
# pour s = 100
## calcul des contributions du pool initial ----
esp_x1 <- param_indiv[1,]$lambda * param_indiv[1,]$alpha / param_indiv[1,]$beta  
esp_x2 <- param_indiv[2,]$lambda * param_indiv[2,]$alpha / param_indiv[2,]$beta 
esp_x3 <- param_indiv[3,]$lambda * param_indiv[3,]$alpha / param_indiv[3,]$beta 

esp_s_init_prorata <- sum(esp_x1, esp_x2, esp_x3)

contrib_x1 <- s * esp_x1 / esp_s_init_prorata
contrib_x2 <- s * esp_x2 / esp_s_init_prorata
contrib_x3 <- s * esp_x3 / esp_s_init_prorata


## calcul des contributions du pool redistribué ----
esp_y1 <- param_indiv[1,]$lambda * param_indiv[1,]$alpha / param_indiv[1,]$beta  
esp_y2 <- sum(fx23 * (seq(nfft) - 1))

esp_s_redist_prorata <- sum(esp_y1, esp_y2)

contrib_y1 <- s * esp_y1 / esp_s_redist_prorata
contrib_y2 <- s * esp_y2 / esp_s_redist_prorata
    
# uniforme ----
# pour s = 1
## calcul des contributions du pool initial ----
pool_size_init <- nrow(param_indiv)

contrib_xi <- s / pool_size_init

## calcul des contributions du pool redistribué ----
pool_size_redist <- nrow(param_indiv) - 1

contrib_yi <- s / pool_size_redist




                   