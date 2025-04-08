###
### Travail 3, ACT-7119
### Exemple numérique : Illustration de la propriété de "volonté à joindre"
###
##


# Importation des paquetages
library(ggplot2)
library(ggpubr)


# Importation des données de sous-pools
subpool1 <- read.csv("donnees-sous-pools/sous-pool-1-k3--2025-04-04.csv")
subpool2 <- read.csv("donnees-sous-pools/sous-pool-2-k3--2025-04-04.csv")
subpool3 <- read.csv("donnees-sous-pools/sous-pool-3-k3--2025-04-04.csv")

source("algorithms/uniform_rsharing.R")
source("algorithms/cond_mean_rsharing_FFT_Blier.R")
source("algorithms/mean_prop_rsharing.R")


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


# Tableau 2
risk_ech_subpool1 <- subpool1_ech[subpool1_ech$id == id_risk_ech_subpool1,]
risk_ech_subpool2 <- subpool2_ech[subpool2_ech$id == id_risk_ech_subpool2,]
risk_ech_subpool3 <- subpool3_ech[subpool3_ech$id == id_risk_ech_subpool3,]


set.seed(11)
subpool123_ech <- rbind(
    risk_ech_subpool1,
    risk_ech_subpool2,
    risk_ech_subpool3,
    subpool1[sample(seq(SUBPOOL1_SIZE), 74),],
    subpool2[sample(seq(SUBPOOL2_SIZE), 74),],
    subpool3[sample(seq(SUBPOOL3_SIZE), 74),]
)

subpool123_ech_SIZE <- nrow(subpool123_ech)


# Échantillon du sous-pool 1
nfft <- 2^17
h <- 1
support_S <- (c(0, seq(nfft - 1)) * h)

ES_pool_melange <- sum(subpool123_ech$lambda * subpool123_ech$alpha / subpool123_ech$beta)
EX_pool_melange_r1 <- subpool123_ech$lambda[1] * subpool123_ech$alpha[1] / subpool123_ech$beta[1]
EX_pool_melange_r2 <- subpool123_ech$lambda[2] * subpool123_ech$alpha[2] / subpool123_ech$beta[2]
EX_pool_melange_r3 <- subpool123_ech$lambda[3] * subpool123_ech$alpha[3] / subpool123_ech$beta[3]

unif_rsharing_pool_melange <- support_S / subpool123_ech_SIZE

mean_prop_rsharing_pool_melange_r1 <- (EX_pool_melange_r1 / ES_pool_melange) * support_S
mean_prop_rsharing_pool_melange_r2 <- (EX_pool_melange_r2 / ES_pool_melange) * support_S
mean_prop_rsharing_pool_melange_r3 <- (EX_pool_melange_r3 / ES_pool_melange) * support_S

cm_rsharing_pool_melange <- compute_cm(nfft, h, subpool123_ech, "unbiased")
sum(cm_rsharing_pool_melange$fs)
# Temps : 56.2 secondes, 2025-04-08, avec "unbiased".
# nfft = 2^17, h = 1.
cm_rsharing_pool_melange_r1 <- cm_rsharing_pool_melange$contrib[[1]]
cm_rsharing_pool_melange_r2 <- cm_rsharing_pool_melange$contrib[[2]]
cm_rsharing_pool_melange_r3 <- cm_rsharing_pool_melange$contrib[[3]]


# Simulation
simulate_pool_losses <- function(n, pool_params)
{
    set.seed(11)
    N_RISKS <- nrow(pool_params)
    realizations <- matrix(numeric(n * N_RISKS), nrow = n)

    # Barre de progression
    progress_bar <- txtProgressBar(1, n, 1, "*", style = 3)
    for (i in seq(n))
    {
        for (j in seq(nrow(pool_params)))
        {
            realizations[i, j] <- simulate_PoisComp_gamma(
                pool_params$lambda[j],
                pool_params$alpha[j],
                pool_params$beta[j]
            )
        }

        setTxtProgressBar(progress_bar, i)
    }
    list(all_losses = realizations,
         total_losses = rowSums(realizations))
}


nsimulations <- 100000
t1 <- Sys.time()
resultats_simulation <- simulate_pool_losses(nsimulations, subpool123_ech)
t2 <- Sys.time()
temps_simulation <- t2 - t1
# 1.87 minutes (n = 100 000 simulations).
# Validation :
c(empirique = mean(resultats_simulation$total_losses),
  theorique = ES_pool_melange) # Ok.

varX1_empirique <- var(resultats_simulation$all_losses[, 1])
varX2_empirique <- var(resultats_simulation$all_losses[, 2])
varX3_empirique <- var(resultats_simulation$all_losses[, 3])


varH1prop <- var((EX_pool_melange_r1 / ES_pool_melange) * resultats_simulation$total_losses)
varH2prop <- var((EX_pool_melange_r2 / ES_pool_melange) * resultats_simulation$total_losses)
varH3prop <- var((EX_pool_melange_r3 / ES_pool_melange) * resultats_simulation$total_losses)
# Ordonnancement dans ce cas, voir [Denuit et Robert, 2021c] dans article de
# base.


varH1uni <- var(resultats_simulation$total_losses / subpool123_ech_SIZE)
varH2uni <- var(resultats_simulation$total_losses / subpool123_ech_SIZE)
varH3uni <- var(resultats_simulation$total_losses / subpool123_ech_SIZE)


ceiling_s_losses <- ceiling(resultats_simulation$total_losses)

varH1cm <- var(cm_rsharing_pool_melange_r1[ceiling_s_losses])
varH2cm <- var(cm_rsharing_pool_melange_r2[ceiling_s_losses])
varH3cm <- var(cm_rsharing_pool_melange_r3[ceiling_s_losses])
