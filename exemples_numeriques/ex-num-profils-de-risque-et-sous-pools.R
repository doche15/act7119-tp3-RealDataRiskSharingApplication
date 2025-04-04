###
### Travail 3, ACT-7119
### Exemple numérique : Profils de risque et sous-pools
###
##


# Importation des paquetages
library(ggplot2)


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


# Tableau 1
summary(subpool1_ech)
summary(subpool2_ech)
summary(subpool3_ech)


# Tableau 2
risk_ech_subpool1 <- subpool1_ech[subpool1_ech$id == id_risk_ech_subpool1,]
risk_ech_subpool2 <- subpool2_ech[subpool2_ech$id == id_risk_ech_subpool2,]
risk_ech_subpool3 <- subpool3_ech[subpool3_ech$id == id_risk_ech_subpool3,]


# Partage de risque
## Restructuration des échantillons de sous-pools pour que le risque
## analysé de chaque sous-pool soit à la ligne 1 des échantillons
subpool1_ech <- rbind(
    risk_ech_subpool1,
    subpool1_ech[!(subpool1_ech$id == id_risk_ech_subpool1),]
)
nrow(subpool1_ech)


subpool2_ech <- rbind(
    risk_ech_subpool2,
    subpool2_ech[!(subpool2_ech$id == id_risk_ech_subpool2),]
)
nrow(subpool2_ech)


subpool3_ech <- rbind(
    risk_ech_subpool3,
    subpool3_ech[!(subpool3_ech$id == id_risk_ech_subpool3),]
)
nrow(subpool3_ech)


## Partage de risque
source("algorithms/uniform_rsharing.R")
source("algorithms/cond_mean_rsharing_FFT_Blier.R")
source("algorithms/mean_prop_rsharing.R")
source("algorithms/linear_reg_rsharing.R")


# Échantillon du sous-pool 1
nfft <- 2^17
h <- 5
support_S <- (c(0, seq(nfft - 1)) * h)

ES_sp1 <- sum(subpool1_ech$lambda * subpool1_ech$alpha / subpool1_ech$beta)
EX_sp1_r1 <- subpool1_ech$lambda[1] * subpool1_ech$alpha[1] / subpool1_ech$beta[1]
VarS_sp1 <- sum((subpool1_ech$lambda * subpool1_ech$alpha / (subpool1_ech$beta^2)) * (subpool1_ech$alpha + 1))
VarX_sp1_r1 <- (subpool1_ech$lambda[1] * subpool1_ech$alpha[1] / (subpool1_ech$beta[1]^2)) * (subpool1_ech$alpha[1] + 1)

unif_rsharing_sp1_r1 <- support_S / RISK_ECH_SIZE
mean_prop_rsharing_sp1_r1 <- (EX_sp1_r1 / ES_sp1) * support_S
linear_reg_rsharing_sp1_r1 <- EX_sp1_r1 + (VarX_sp1_r1 / VarS_sp1) * (support_S - ES_sp1)
cm_rsharing_sp1 <- compute_cm(nfft, h, subpool1_ech, "unbiased")
# Temps : 249.3 secondes, 2025-04-04, avec "unbiased".
cm_rsharing_sp1_r1 <- cm_rsharing_sp1$contrib[[1]]

results_contrib_rsharing_sp1 <- data.frame(s = support_S,
                                           unif_rsharing = unif_rsharing_sp1_r1,
                                           mean_prop_rsharing = mean_prop_rsharing_sp1_r1,
                                           linear_reg_rsharing = linear_reg_rsharing_sp1_r1,
                                           cm_rsharing = cm_rsharing_sp1_r1)

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, cm_rsharing_sp1$fs)
sum(cm_rsharing_sp1$fs)
# Entre environ 30 000 et 90 000.

cm_rsharing_sp1_r1[30100 / h]
cm_rsharing_sp1_r1[89900 / h]

results_contrib_rsharing_sp1_reduced <- results_contrib_rsharing_sp1[(30100:89900) / h,]

results_contrib_rsharing_sp1_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp1_reduced$s,
                                                                contribution = results_contrib_rsharing_sp1_reduced$unif_rsharing,
                                                                methode = rep("Uniforme", nrow(results_contrib_rsharing_sp1_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp1_reduced$s,
                                                                contribution = results_contrib_rsharing_sp1_reduced$mean_prop_rsharing,
                                                                methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp1_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp1_reduced$s,
                                                                contribution = results_contrib_rsharing_sp1_reduced$linear_reg_rsharing,
                                                                methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp1_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp1_reduced$s,
                                                                contribution = results_contrib_rsharing_sp1_reduced$cm_rsharing,
                                                                methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp1_reduced))))


# Référence :
# https://stackoverflow.com/questions/14794599/how-to-change-line-width-in-ggplot
ggplot(results_contrib_rsharing_sp1_reduced_joined) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    xlab("s (en milliers)") +
    ylab("Contribution") +
    xlim(30000 / 1000, 90000 / 1000) +
    ylim(20, 100) +
    theme_classic() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          legend.text = element_text(size=9),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))


# Échantillon du sous-pool 2
ES_sp2 <- sum(subpool2_ech$lambda * subpool2_ech$alpha / subpool2_ech$beta)
EX_sp2_r1 <- subpool2_ech$lambda[1] * subpool2_ech$alpha[1] / subpool2_ech$beta[1]
VarS_sp2 <- sum((subpool2_ech$lambda * subpool2_ech$alpha / (subpool2_ech$beta^2)) * (subpool2_ech$alpha + 1))
VarX_sp2_r1 <- (subpool2_ech$lambda[1] * subpool2_ech$alpha[1] / (subpool2_ech$beta[1]^2)) * (subpool2_ech$alpha[1] + 1)

unif_rsharing_sp2_r1 <- support_S / RISK_ECH_SIZE
mean_prop_rsharing_sp2_r1 <- (EX_sp2_r1 / ES_sp2) * support_S
linear_reg_rsharing_sp2_r1 <- EX_sp2_r1 + (VarX_sp2_r1 / VarS_sp2) * (support_S - ES_sp2)
cm_rsharing_sp2 <- compute_cm(nfft, h, subpool2_ech, "unbiased")
# Temps : 246.3 secondes, 2025-04-04, avec "unbiased".
# Erreur obtenue : "Error in order(typeScores, scores) : argument lengths differ"
cm_rsharing_sp2_r1 <- cm_rsharing_sp2$contrib[[1]]

results_contrib_rsharing_sp2 <- data.frame(s = support_S,
                                           unif_rsharing = unif_rsharing_sp2_r1,
                                           mean_prop_rsharing = mean_prop_rsharing_sp2_r1,
                                           linear_reg_rsharing = linear_reg_rsharing_sp2_r1,
                                           cm_rsharing = cm_rsharing_sp2_r1)

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, cm_rsharing_sp2$fs)
sum(cm_rsharing_sp2$fs)
# Entre environ 50 000 et 120 000.

cm_rsharing_sp2_r1[50000 / h]
cm_rsharing_sp2_r1[120000 / h]


results_contrib_rsharing_sp2_reduced <- results_contrib_rsharing_sp2[(50100:119900) / h,]

results_contrib_rsharing_sp2_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                                                                contribution = results_contrib_rsharing_sp2_reduced$unif_rsharing,
                                                                methode = rep("Uniforme", nrow(results_contrib_rsharing_sp2_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                                                                contribution = results_contrib_rsharing_sp2_reduced$mean_prop_rsharing,
                                                                methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp2_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                                                                contribution = results_contrib_rsharing_sp2_reduced$linear_reg_rsharing,
                                                                methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp2_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                                                                contribution = results_contrib_rsharing_sp2_reduced$cm_rsharing,
                                                                methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp2_reduced))))


# Référence :
# https://stackoverflow.com/questions/14794599/how-to-change-line-width-in-ggplot
ggplot(results_contrib_rsharing_sp2_reduced_joined) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    xlab("s (en milliers)") +
    ylab("Contribution") +
    xlim(50000 / 1000, 120000 / 1000) +
    ylim(50, 130) +
    theme_classic() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          legend.text = element_text(size=9),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))


# Échantillon du sous-pool 3
ES_sp3 <- sum(subpool3_ech$lambda * subpool3_ech$alpha / subpool3_ech$beta)
EX_sp3_r1 <- subpool3_ech$lambda[1] * subpool3_ech$alpha[1] / subpool3_ech$beta[1]
VarS_sp3 <- sum((subpool3_ech$lambda * subpool3_ech$alpha / (subpool3_ech$beta^2)) * (subpool3_ech$alpha + 1))
VarX_sp3_r1 <- (subpool3_ech$lambda[1] * subpool3_ech$alpha[1] / (subpool3_ech$beta[1]^2)) * (subpool3_ech$alpha[1] + 1)

unif_rsharing_sp3_r1 <- support_S / RISK_ECH_SIZE
mean_prop_rsharing_sp3_r1 <- (EX_sp3_r1 / ES_sp3) * support_S
linear_reg_rsharing_sp3_r1 <- EX_sp3_r1 + (VarX_sp3_r1 / VarS_sp3) * (support_S - ES_sp3)
cm_rsharing_sp3 <- compute_cm(nfft, h, subpool3_ech, "unbiased")
# Temps : 258.1 secondes, 2025-04-04, avec "unbiased".
cm_rsharing_sp3_r1 <- cm_rsharing_sp3$contrib[[1]]

results_contrib_rsharing_sp3 <- data.frame(s = support_S,
                                           unif_rsharing = unif_rsharing_sp3_r1,
                                           mean_prop_rsharing = mean_prop_rsharing_sp3_r1,
                                           linear_reg_rsharing = linear_reg_rsharing_sp3_r1,
                                           cm_rsharing = cm_rsharing_sp3_r1)

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, cm_rsharing_sp3$fs)
sum(cm_rsharing_sp3$fs)
# Entre environ 110 000 et 220 000.

cm_rsharing_sp3_r1[110000 / h]
cm_rsharing_sp3_r1[220000 / h]


results_contrib_rsharing_sp3_reduced <- results_contrib_rsharing_sp3[(110100:219900) / h,]

results_contrib_rsharing_sp3_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp3_reduced$s,
                                                                contribution = results_contrib_rsharing_sp3_reduced$unif_rsharing,
                                                                methode = rep("Uniforme", nrow(results_contrib_rsharing_sp3_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp3_reduced$s,
                                                                contribution = results_contrib_rsharing_sp3_reduced$mean_prop_rsharing,
                                                                methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp3_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp3_reduced$s,
                                                                contribution = results_contrib_rsharing_sp3_reduced$linear_reg_rsharing,
                                                                methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp3_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp3_reduced$s,
                                                                contribution = results_contrib_rsharing_sp3_reduced$cm_rsharing,
                                                                methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp3_reduced))))


# Référence :
# https://stackoverflow.com/questions/14794599/how-to-change-line-width-in-ggplot
ggplot(results_contrib_rsharing_sp3_reduced_joined) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    xlab("s (en milliers)") +
    ylab("Contribution") +
    xlim(110000 / 1000, 220000 / 1000) +
    ylim(50, 300) +
    theme_classic() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          legend.text = element_text(size=9),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))




# Validation
# results_contrib_rsharing_sp3_reduced <- results_contrib_rsharing_sp3[(80000:250000) / h,]
#
# results_contrib_rsharing_sp3_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$unif_rsharing,
#                                                                 methode = rep("Uniforme", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$mean_prop_rsharing,
#                                                                 methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$linear_reg_rsharing,
#                                                                 methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$cm_rsharing,
#                                                                 methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp3_reduced))))
# ggplot(results_contrib_rsharing_sp3_reduced_joined) +
#     geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
#     xlab("s (en milliers)") +
#     ylab("Contribution") +
#     xlim(80000 / 1000, 250000 / 1000) +
#     ylim(50, 300) +
#     theme_classic() +
#     theme(axis.text = element_text(size=12),
#           axis.title = element_text(size=12),
#           legend.text = element_text(size=9),
#           legend.position = "bottom") +
#     guides(col = guide_legend(title = "Méthode"))
